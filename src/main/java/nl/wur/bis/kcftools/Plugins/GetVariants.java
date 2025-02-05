package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.atomic.AtomicInteger;

@Command(name = "getVariations", description = " Screen for reference kmers that are not present in the KMC database, and detect variation")
public class GetVariants implements Callable<Integer>, Runnable {
    // in reference file name
    @Option(names = {"-r", "--reference"}, description = "Reference file name", required = true)
    private String refFasta;
    // in KMC database prefix
    @Option(names = {"-k", "--kmc"}, description = "KMC database prefix", required = true)
    private String kmcDBprefix;
    // in window size
    @Option(names = {"-w", "--window"}, description = "Window size", required = false)
    private int windowSize;
    // in output file name
    @Option(names = {"-o", "--output"}, description = "Output file name", required = true)
    private String outFile;
    // in number of threads
    @Option(names = {"-t", "--threads"}, description = "Number of threads", required = false)
    private int nThreads = 2;
    // sample name
    @Option(names = {"-s", "--sample"}, description = "Sample name", required = true)
    private String sampleName;
    // load kmc into memory
    @Option(names = {"-m", "--memory"}, description = "Load KMC database into memory", required = false)
    private boolean loadMemory = false;
    // if its RNA-seq data, get gtf file and feature type (gene, mRNA, CDS)
    @Option(names = {"-g", "--gtf"}, description = "GTF file name", required = false)
    private String gtfFile;
    // feature type: gene or transcript
    @Option(names = {"-f", "--feature"}, description = "Feature type (\"gene\" or \"transcript\")", required = false)
    private String featureType;

    private final String CLASS_NAME = this.getClass().getSimpleName();
    private String model = "wholegenome";
    private FastaIndex index;
    private int kmerSize;
    private GTFReader gtfReader;

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
     * Get the variations
     */
    public void getVariations() throws IOException {

        KMC kmc = new KMC(kmcDBprefix, loadMemory);
        kmerSize = kmc.getKmerLength();
        KCFHeader header = new KCFHeader();
        header.setReference(refFasta);
        header.addCommandLine(HelperFunctions.getCommandLine());
        header.addSample(sampleName);
        header.setWindowSize(windowSize);
        header.setKmerSize(kmc.getKmerLength());
        header.setIBS(false);

        index = new FastaIndex(refFasta);
        ConcurrentHashMap<String, Queue<Window>> windowsMap = new ConcurrentHashMap<>();

        if (gtfFile != null) {
            model = "targeted";
            gtfReader = new GTFReader(gtfFile);
        }
        else{
            model = "wholegenome";
        }

        HelperFunctions.log("info", CLASS_NAME, "Generating windows...");
        for (String name : index.getSequenceNames()) {
            header.addContig(name, index.getSequenceLength(name));
            Queue<Window> windows = getWindows(name);
            windowsMap.put(name, windows);
        }

        int totalWindows = windowsMap.values().stream().mapToInt(Queue::size).sum();
        HelperFunctions.log("info", CLASS_NAME, "Number of windows: " + totalWindows);

        LinkedHashMap<String, List<Window>> processedWindows = new LinkedHashMap<>();
        AtomicInteger completedWindows = new AtomicInteger(0);

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        for (String name : windowsMap.keySet()) {
            Queue<Window> windows = windowsMap.get(name);
            List<Window> processed = new ArrayList<>(windows.size());

            ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);
            for (Window window : windows) {
                completionService.submit(() -> {
                    Fasta fasta = getFasta(window);
                    Window processedWindow = processWindow(window, fasta, kmc);
                    synchronized (processed) {
                        processed.add(processedWindow);
                    }
                    int completed = completedWindows.incrementAndGet();
                    float progress = (float) (completed * 100) / totalWindows;
                    synchronized (System.out) {
                        System.out.printf("\rProgress: %.2f%%", progress);
                    }
                    return null;
                });
            }

            // wait for all the windows to be processed
            for (int i = 0; i < windows.size(); i++) {
                try {
                    // let's wait to get the results
                    completionService.take().get();
                } catch (InterruptedException | ExecutionException e) {
                    e.printStackTrace();
                }
            }
            processedWindows.put(name, processed);
        }
        executor.shutdown();

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
        return switch (model) {
            case "wholegenome" -> window.getFasta(index);
            case "targeted" -> gtfReader.getFasta(window.getWindowId(), featureType, index);
            default -> {
                HelperFunctions.log("error", CLASS_NAME, "Invalid model type: " + model + ". Supported models are 'wholegenome' and 'targeted'");
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
//        int localDistance = 0;
        int localInnerDistance = 0;
        int localTailDistance = 0;
        int gapSize = 0;
        boolean isTail = true;
//        int tailGap = 0;

        if (fasta == null) {
            HelperFunctions.log("error", CLASS_NAME, "Fasta object is null for window: " + window.getWindowId());
            return window;
        }
        List<Kmer> kmers = fasta.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);

        if (!kmers.isEmpty()) {
            for (Kmer k : kmers) {
                localTotalKmers++;
                Kmer km = new Kmer(k, kmc.isBothStrands());

                if (kmc.isExist(km)) {
                    // if the kmer exists in the KMC database, 1+ observed kmers
                    localObservedKmers++;
                    if (gapSize > 0) {
                        // if there is a gap, increment the variation
                        localVariation++;
                        // if the gap is at the beginning or end of the window, increment the distance
                        if (isTail) {
                            localTailDistance += gapSize;
//                            localDistance += gapSize;
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

            // Process the last gap if it exists
            if (gapSize > 0) {
                localVariation++;
//                localDistance += gapSize;
                localTailDistance += gapSize;
            }
        }

        synchronized (window) {
            window.addTotalKmers(localTotalKmers);
            window.setEffLength(fasta.getEffectiveATGCCount(kmerSize));
            window.addData(sampleName, localObservedKmers, localVariation, localInnerDistance, localTailDistance, "N");
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
        switch (model) {
            case "wholegenome" -> {
                int lastEnd = 0;
                while (lastEnd < sequenceLength) {
                    int start = Math.max(0, lastEnd - kmerSize);
                    int end = Math.min(start + windowSize, sequenceLength);
                    if (end - start >= kmerSize) {
                        windows.add(new Window(sequenceName + "_" + start, sequenceName, start, end));
                    }
                    lastEnd = end;
                }
            }
            case "targeted" -> {
                for (GTFReader.Feature feature : gtfReader.getFeatures(sequenceName, featureType)) {
                    windows.add(new Window(feature.getId(), sequenceName, feature.getStart(), feature.getEnd()));
                }
            }
            default -> HelperFunctions.log("error", CLASS_NAME, "Invalid model type: " + model + ". Supported models are 'wholegenome' and 'targeted'");
        }
        return windows;
    }

    /***
     * Validate the command line arguments
     */
    private void validateCMD() {
        if(gtfFile != null && gtfFile.isEmpty()){
            if (windowSize > 0) {
                HelperFunctions.log("error", CLASS_NAME, "Window size is not valid for targeted model");
            }
            if (featureType == null) {
                HelperFunctions.log("error", CLASS_NAME, "Feature type (gene or transcript) is required for targeted model");
            }
            if (!featureType.equals("gene") && !featureType.equals("transcript")) {
                HelperFunctions.log("error", CLASS_NAME, "Invalid feature type: " + featureType + ". Supported feature types are 'gene' and 'transcript'");
            }
        }
        if (featureType != null && !featureType.isEmpty()) {
            if (gtfFile == null) {
                HelperFunctions.log("error", CLASS_NAME, "GTF file is required for targeted model");
            }
        }
        if (windowSize <= 0 && gtfFile == null) {
            HelperFunctions.log("error", CLASS_NAME, "Window size (for wholegenome model) or GTF file (for targeted model) is required");
        }
        if (nThreads <= 0) {
            HelperFunctions.log("error", CLASS_NAME, "Number of threads should be greater than 0");
        }
    }
}
//EOF