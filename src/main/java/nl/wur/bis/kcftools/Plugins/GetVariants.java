package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;

import picocli.CommandLine;
import picocli.CommandLine.*;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

/***
 * This is a command line plugin that may extract kmers from the reference file and compare them to the KMC database
 * It will detect the variations in the kmers and write them to a KCF file
 */
@Command(name = "getVariations", description = " Screen for reference kmers that are not present in the KMC database, and detect variation", sortOptions = false)
public class GetVariants implements Callable<Integer>, Runnable {
    // in reference file name
    @Option(names = {"-r", "--reference"}, description = "Reference file name", required = true)
    private String refFasta;
    // in KMC database prefix
    @Option(names = {"-k", "--kmc"}, description = "KMC database prefix", required = true)
    private String kmcDBprefix;
    // in output file name
    @Option(names = {"-o", "--output"}, description = "Output file name", required = true)
    private String outFile;
    // sample name
    @Option(names = {"-s", "--sample"}, description = "Sample name", required = true)
    private String sampleName;
    // feature type: window or gene or transcript
    @Option(names = {"-f", "--feature"}, description = "Feature type (\"window\" or \"gene\" or \"transcript\")", required = true)
    private String featureType;
    // in number of threads
    @Option(names = {"-t", "--threads"}, description = "Number of threads [2]", required = false)
    private int nThreads = 2;
    // load kmc into memory
    @Option(names = {"-m", "--memory"}, description = "Load KMC database into memory", required = false)
    private boolean loadMemory = false;
    // inner kmer distance weight
    @Option(names = {"--wi"}, description = "Inner kmer distance weight [0.3]", required = false)
    private double innerDistanceWeight = 0.3;
    // tail kmer distance weight
    @Option(names = {"--wt"}, description = "Tail kmer distance weight [0.3]", required = false)
    private double tailDistanceWeight = 0.3;
    // kmer ratio weight
    @Option(names = {"--wr"}, description = "Kmer ratio weight [0.4]", required = false)
    private double kmerRatioWeight = 0.4;
    // in window size
    @Option(names = {"-w", "--window"}, description = "Window size", required = false)
    private int windowSize;
    // if its RNA-seq data, get gtf file and feature type (gene, mRNA, CDS)
    @Option(names = {"-g", "--gtf"}, description = "GTF file name", required = false)
    private String gtfFile;
    // minimum kmer count to consider a kmer as valid
    @Option(names = {"-c", "--min-k-count"}, description = "Minimum kmer count to consider [1]", required = false)
    private int minKmerCount = 1;
    // step size for sliding window
    @Option(names = {"-p", "--step"}, description = "Step size for sliding window [window size]", required = false)
    private int stepSize = 0;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private FastaIndex index;
    private int kmerSize;
    private GTF gtf;
//    private final double[] weights = new double[] {innerDistanceWeight, tailDistanceWeight, kmerRatioWeight};

    public GetVariants() {
    }

    @Override
    public Integer call() throws IOException, InterruptedException {
        HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        validateCMD();
        getVariations();
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (IOException | InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    /***
     * Main function to get the variations
     */
    public void getVariations() throws IOException {

        sampleName = cleanSampleName(sampleName);
        KMC kmc = new KMC(kmcDBprefix, loadMemory);
        kmerSize = kmc.getKmerLength();
        KCFHeader header = new KCFHeader();
        header.setReference(refFasta);
        header.addCommandLine(HelperFunctions.getCommandLine());
        header.addSample(sampleName);
        header.setWindowSize(windowSize);
        header.setStepSize(stepSize);
        header.setKmerSize(kmc.getKmerLength());
        header.setIBS(false);
        header.setWeightInnerDist(innerDistanceWeight);
        header.setWeightTailDist(tailDistanceWeight);
        header.setWeightKmerRatio(kmerRatioWeight);

        index = new FastaIndex(refFasta);
        ConcurrentHashMap<String, Queue<Window>> windowsMap = new ConcurrentHashMap<>();

        if (featureType.equals("gene") || featureType.equals("transcript")){
            gtf = new GTF(gtfFile);
        }

        Logger.info(CLASS_NAME, "Generating windows...");
        for (String name : index.getSequenceNames()) {
            header.addContig(name, index.getSequenceLength(name));
            Queue<Window> windows = getWindows(name);
            windowsMap.put(name, windows);
        }

        int totalWindows = windowsMap.values().stream().mapToInt(Queue::size).sum();
        Logger.info(CLASS_NAME, "Number of windows: " + totalWindows);

        LinkedHashMap<String, List<Window>> processedWindows = new LinkedHashMap<>();
        AtomicInteger completedWindows = new AtomicInteger(0);

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);

        for (Map.Entry<String, Queue<Window>> entry : windowsMap.entrySet()) {
            String name = entry.getKey();
            Queue<Window> windows = entry.getValue();
            List<Window> processed = Collections.synchronizedList(new ArrayList<>(windows.size()));
            processedWindows.put(name, processed);

            for (Window window : windows) {
                completionService.submit(() -> {
                    Fasta fasta = getFasta(window);
                    Window processedWindow = processWindow(window, fasta, kmc);
                    processed.add(processedWindow);
                    int completed = completedWindows.incrementAndGet();
                    float progress = (float) (completed * 100) / totalWindows;
                    synchronized (System.out) {
                        System.out.printf("\rProgress: %.2f%%", progress);
                    }
                    return null;
                });
            }
        }

        for (int i = 0; i < totalWindows; i++) {
            try {
                completionService.take().get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        System.out.print("\r");
        for (int i = 0; i < 100; i++) {
            System.out.print(" ");
        }
        System.out.print("\r");
        System.out.flush();

        // sort the windows based on its start position
        for (String name : processedWindows.keySet()) {
            processedWindows.get(name).sort(Comparator.comparingInt(Window::getStart));
        }

        header.setWindowCount(totalWindows);
        try (KCFWriter writer = new KCFWriter(outFile)) {
            writer.writeHeader(header);
            for (String name : index.getSequenceNames()) {
                writer.writeWindows(processedWindows.get(name));
            }
        }
        index.close();
        kmc.close();
        HelperFunctions.printMaxMemoryUsage();
    }

    /***
     * Get the Fasta object based on the model type
     */
    private Fasta getFasta(Window window){
        return switch (featureType) {
            case "window" -> window.getFasta(index);
            case "gene", "transcript" -> gtf.getFasta(window.getWindowId(), index, featureType.equals("gene"));
            default -> {
                Logger.error(CLASS_NAME, "Invalid model type: " + featureType + ". Supported models are 'window' or 'gene' or 'transcript'");
                yield null;
            }
        };
    }

    /***
     * Process a window and calculate the number of observed kmers and the variation
     */
    private Window processWindow(Window window, Fasta fasta, KMC kmc) {
        int localTotalKmers = 0;
        int localObservedKmers = 0;
        int localVariation = 0;
        int localInnerDistance = 0;
        int gapSize = 0;
        boolean isTail = true;
        int localLeftDist = 0;
        int localRightDist = 0;
        long localKmerCount = 0;

        if (fasta == null) {
            Logger.error(CLASS_NAME, "Fasta object is null for window: " + window.getWindowId());
            return window;
        }
        List<Kmer> kmers = fasta.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);

        if (!kmers.isEmpty()) {
            for (Kmer k : kmers) {
                localTotalKmers++;
                Kmer km = new Kmer(k, kmc.isBothStrands());
                int kmerCount = kmc.getCount(km);
                if (kmerCount >= minKmerCount) {
                    // if the kmer exists in the KMC database, 1+ observed kmers
                    localKmerCount += kmerCount;
                    localObservedKmers++;
                    if (gapSize > 0) {
                        // if there is a gap, increment the variation
                        localVariation++;
                        // if the gap is at the beginning or end of the window, increment the distance
                        if (isTail) {
                            localLeftDist += gapSize;
                        }
                        else {
                            // if the gap is in the middle of the window, calculate the distance based on the gap size and kmer size
                            localInnerDistance += getDistance(kmc, gapSize);
                        }
                    }
                    isTail = false;
                    gapSize = 0;
                } else {
                    gapSize++;
                }
            }

            // process the last gap if it exists
            if (gapSize > 0) {
                localVariation++;
                localRightDist += gapSize;
            }
        }

        synchronized (window) {
            window.addTotalKmers(localTotalKmers);
            window.setEffLength(fasta.getEffectiveATGCCount(kmerSize));
            window.addData(sampleName, localObservedKmers, localVariation, localInnerDistance, localLeftDist, localRightDist, localKmerCount, "N", getWeights());
        }

        return window;
    }

    /***
     * Get the distance based on the gap size
     * if the actual missing base is 1, and the kmer size is 3, we will have 3 missing kmers, hence the distance is 1 (3 - (3-1)) = 1
     */
    private static int getDistance(KMC kmc, int gapSize) {
        int distance = gapSize - (kmc.getKmerLength() - 1);
        if (distance <= 0) {
            distance = Math.abs(distance + 1);
        }
        return distance;
    }

    /***
     * Get the windows based on the model type
     */
    private Queue<Window> getWindows(String sequenceName) {
        Queue<Window> windows = new ConcurrentLinkedDeque<>();
        int sequenceLength = index.getSequenceLength(sequenceName);
        switch (featureType) {
            case "window" -> {
//                int lastEnd = 0;
//                while (lastEnd < sequenceLength) {
//                    int start = Math.max(0, lastEnd - kmerSize + 1);
//                    int end = Math.min(start + windowSize, sequenceLength);
//                    if (end - start >= kmerSize) {
//                        windows.add(new Window(sequenceName + "_" + start, sequenceName, start, end));
//                    }
//                    lastEnd = end;
//                }

                int lastPos = 0;

                if (stepSize > 0) {
                    // Sliding window mode
                    while (lastPos < sequenceLength) {
                        int start = lastPos;
                        int end = Math.min(start + windowSize, sequenceLength);

                        if (end - start >= kmerSize) {
                            windows.add(new Window(sequenceName + "_" + start, sequenceName, start, end));
                        }

                        lastPos += stepSize; // slide by stepSize
                    }
                } else {
                    // Old tiling mode (non-overlapping but with kmer overlap at boundary)
                    int lastEnd = 0;
                    while (lastEnd < sequenceLength) {
                        int start = Math.max(0, lastEnd - kmerSize + 1);
                        int end = Math.min(start + windowSize, sequenceLength);

                        if (end - start >= kmerSize) {
                            windows.add(new Window(sequenceName + "_" + start, sequenceName, start, end));
                        }

                        lastEnd = end;
                    }
                }


            }
            case "gene" -> {
                String[] genes = gtf.getGenes(sequenceName);
                for (String gene : genes) {
                    GTF.Loci loci = gtf.getLoci(gene);
                    windows.add(new Window(gene, loci.getChromosome(), loci.getStart(), loci.getEnd()));
                }
            }
            case "transcript" -> {
                String[] genes = gtf.getGenes(sequenceName);
                if (genes.length == 0) {
                    Logger.warning(CLASS_NAME, "No genes found in GTF file for sequence: " + sequenceName);
                    return windows;
                }
                for (String gene : genes) {
                    String[] transcripts = gtf.getTranscripts(gene);
                    if (transcripts.length == 0) {
                        Logger.error(CLASS_NAME, "No transcripts found for gene: " + gene + " in GTF file for sequence: " + sequenceName);
                        continue;
                    }
                    for (String transcript : transcripts) {
                        GTF.Loci loci = gtf.getLoci(transcript);
                        windows.add(new Window(transcript, loci.getChromosome(), loci.getStart(), loci.getEnd()));
                    }
                }
            }
            default -> Logger.error(CLASS_NAME, "Invalid model type: " + featureType + ". Supported models are 'window' or 'gene' or 'transcript'");
        }
        return windows;
    }

    /***
     * Validate the command line arguments
     */
    private void validateCMD() {
        switch (featureType) {
            case "window" -> {
                if (windowSize <= 0) {
                    Logger.error(CLASS_NAME, "Window size is required for window model");
                }
                if (gtfFile != null && !gtfFile.isEmpty()) {
                    Logger.error(CLASS_NAME, "GTF file is not valid for window model");
                }
            }
            case "gene", "transcript" -> {
                if (gtfFile == null || gtfFile.isEmpty()) {
                    Logger.error(CLASS_NAME, "GTF file is required for targeted model");
                }
                if (windowSize > 0) {
                    Logger.error(CLASS_NAME, "Window size is not valid for targeted model");
                }
            }
            default ->
                    Logger.error(CLASS_NAME, "Invalid model type: " + featureType + ". Supported models are 'window' or 'gene' or 'transcript'");
        }

        if (nThreads <= 0) {
            Logger.error(CLASS_NAME, "Number of threads should be greater than 0");
        }

        if (minKmerCount < 1) {
            Logger.error(CLASS_NAME, "Minimum kmer count should be at least 1");
        }
    }

    private double[] getWeights(){
        return new double[] {innerDistanceWeight, tailDistanceWeight, kmerRatioWeight};
    }

    private String cleanSampleName(String sampleName) {
        // Replace invalid filename characters with '_'
        String sanitized = sampleName.replaceAll("[\\\\/:*?\"<>|]", "_");

        if (!sanitized.equals(sampleName)) {
            Logger.warning(CLASS_NAME,
                    "Sample name contains invalid characters, changed to: " + sanitized);
        }
        return sanitized;
    }
}
//EOF