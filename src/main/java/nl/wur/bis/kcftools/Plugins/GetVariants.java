package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
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
    @Option(names = {"-w", "--window"}, description = "Window size", required = true)
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

    private final String CLASS_NAME = this.getClass().getSimpleName();


    public GetVariants() {
    }

    public GetVariants(String refFasta, String kmcDBprefix, int windowSize, String outFile, int nThreads) {
        this.refFasta = refFasta;
        this.kmcDBprefix = kmcDBprefix;
        this.windowSize = windowSize;
        this.outFile = outFile;
        this.nThreads = nThreads;
    }

    @Override
    public Integer call() throws IOException, InterruptedException {
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

    public void getVariations() throws IOException, InterruptedException {
        // eg cmd line: kcftools get_variations -r /data/arabidopsis.fna -k /data/arabidopsis -w 1000 -o /data/arabidopsis.variations.tsv -t 4
        printCommandLine();
        KMC kmc = new KMC(kmcDBprefix, loadMemory);
        KCFHeader header = new KCFHeader();
        header.setReference(refFasta);
        header.addCommandLine(HelperFunctions.getCommandLine());
        header.addSample(sampleName);
        header.setWindowSize(windowSize);
        header.setKmerSize(kmc.getKmerLength());
        header.setIBS(false);

        FastaIndex index = new FastaIndex(refFasta);
        AtomicInteger windowId = new AtomicInteger(0); // For unique window IDs
        ConcurrentHashMap<String, Queue<Window>> windowsMap = new ConcurrentHashMap<>();

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        List<Future<Void>> futures = new ArrayList<>();

        HelperFunctions.log("info", CLASS_NAME, "Generating windows...");
        for (String name : index.getSequenceNames()) {
            header.addContig(name, index.getSequenceLength(name));
            Queue<Window> windows = getWindows(name, index.getSequenceLength(name), kmc.getKmerLength(), windowSize);
            windowsMap.put(name, windows);
        }

        HelperFunctions.log("info", CLASS_NAME, "Update Window ids...");
        // change the windowID for each window in increasing order
        for (String name : index.getSequenceNames()) {
            Queue<Window> windows = windowsMap.get(name);
            for (Window window : windows) {
                window.setWindowId(windowId.incrementAndGet());
            }
        }

        int totalWindows = windowsMap.values().stream().mapToInt(Queue::size).sum();
        HelperFunctions.log("info", CLASS_NAME, "Number of windows: " + totalWindows);

        LinkedHashMap<String, List<Window>> processedWindows = new LinkedHashMap<>();
        AtomicInteger completedWindows = new AtomicInteger(0);

        executor = Executors.newFixedThreadPool(nThreads);
        for (String name : windowsMap.keySet()) {
            Queue<Window> windows = windowsMap.get(name);
            List<Window> processed = new ArrayList<>(windows.size());

            ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);
            for (Window window : windows) {
                completionService.submit(() -> {
                    Window processedWindow = processWindow(window, index, kmc);
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
        printMaxMemoryUsage();
    }




    /***
     * Process a window and calculate the number of observed kmers and the variation
     */
//    private Window processWindow(Window window, FastaIndex index, KMC kmc) {
//        Fasta fasta = window.getFasta(index);
//        int localTotalKmers = 0;
//        int localObservedKmers = 0;
//        int localVariation = 0;
//        int localDistance = 0;
//        int gapSize = 0;
//
//        List<Kmer> kmers = fasta.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);
//        if (!kmers.isEmpty()) {
//            for (Kmer k : kmers) {
//                localTotalKmers++;
//                Kmer km = new Kmer(k, kmc.isBothStrands());
//                if (kmc.isExist(km)) {
//                    localObservedKmers++;
//                    if (gapSize > 0) {
//                        localVariation++;
//                        int distance = gapSize - kmc.getKmerLength() - 1;
//                        if (distance < 0) {
//                            distance = (distance + 1) * -1;
//                        }
//                        localDistance += distance;
//                    }
//                    gapSize = 0;
//                } else {
//                    gapSize++;
//                }
//            }
//            if (gapSize > 0) {
//                localVariation++;
//                int distance = gapSize - (kmc.getKmerLength() - 1);
//                if (distance <= 0) {
//                    distance = Math.abs(distance + 1);
//                }
//                localDistance += distance;
//                gapSize = 0;
//            }
//        }
//        synchronized (window) {
//            window.addTotalKmers(localTotalKmers);
//            window.addData(sampleName, localObservedKmers, localVariation, localDistance, "N");
//        }
//        return window;
//    }

    private Window processWindow(Window window, FastaIndex index, KMC kmc) {
        Fasta fasta = window.getFasta(index);
        int localTotalKmers = 0;
        int localObservedKmers = 0;
        int localVariation = 0;
        int localDistance = 0;
        int gapSize = 0;

        List<Kmer> kmers = fasta.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);
        if (!kmers.isEmpty()) {
            for (Kmer k : kmers) {
                localTotalKmers++;
                Kmer km = new Kmer(k, kmc.isBothStrands());

                if (kmc.isExist(km)) {
                    localObservedKmers++;
                    if (gapSize > 0) {
                        localVariation++;
                        int distance = gapSize - (kmc.getKmerLength() - 1);
                        if (distance <= 0) {
                            distance = Math.abs(distance + 1);
                        }
                        localDistance += distance;
                    }
                    gapSize = 0;
                } else {
                    gapSize++;
                }
            }

            // Process the last gap if it exists
            if (gapSize > 0) {
                localVariation++;
                int distance = gapSize - (kmc.getKmerLength() - 1);
                if (distance <= 0) {
                    distance = Math.abs(distance + 1);
                }
                localDistance += distance;
            }
        }

        synchronized (window) {
            window.addTotalKmers(localTotalKmers);
            window.addData(sampleName, localObservedKmers, localVariation, localDistance, "N");
        }

        return window;
    }


    /***
     * Print the command line options
     */
    private void printCommandLine() {
        HelperFunctions.log("info", CLASS_NAME, "========== CMD options - getVariations ===========");
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %s", "Reference file name", refFasta));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %s", "KMC database prefix", kmcDBprefix));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %d", "Window size", windowSize));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %s", "Output file name", outFile));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %d", "Number of threads", nThreads));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %s", "Sample name", sampleName));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %s", "Load KMC into memory", loadMemory));
        HelperFunctions.log("info", CLASS_NAME, "==================================================");
    }

    private Queue<Window> getWindows(String sequenceName, int sequenceLength, int kmerSize, int windowSize) {
        Queue<Window> windows = new ConcurrentLinkedDeque<>();
        int lastEnd = 0;
        while (lastEnd < sequenceLength) {
            int start = Math.max(0, lastEnd - kmerSize);
            int end = Math.min(start + windowSize, sequenceLength);
            if (end - start >= kmerSize) {
                windows.add(new Window(1, sequenceName, start, end));
            }
            lastEnd = end;
        }
        return windows;
    }

    private void printMaxMemoryUsage() {
        Runtime runtime = Runtime.getRuntime();
        long maxMemory = runtime.maxMemory();
        long allocatedMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();
        long usedMemory = allocatedMemory - freeMemory;

        HelperFunctions.log("info", CLASS_NAME, "============= Memory Usage Statistics ============");
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %.2f", "Max Memory (GB)", maxMemory / (1024.0 * 1024 * 1024)));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %.2f", "Allocated Memory (GB)", allocatedMemory / (1024.0 * 1024 * 1024)));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %.2f", "Free Memory (GB)", freeMemory / (1024.0 * 1024 * 1024)));
        HelperFunctions.log("info", CLASS_NAME, String.format("%-25s: %.2f", "Used Memory (GB)", usedMemory / (1024.0 * 1024 * 1024)));
        HelperFunctions.log("info", CLASS_NAME, "==================================================");
    }
}
//EOF