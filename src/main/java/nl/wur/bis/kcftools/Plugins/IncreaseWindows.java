package nl.wur.bis.kcftools.Plugins;

import java.util.*;
import java.util.concurrent.*;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

/***
 * This is a command line plugin that increases the window size of a KCF file by merging windows
 */
@Command(name = "increaseWindow", description = "Increase the window size of a KCF file by merging windows")
public class IncreaseWindows implements Runnable, Callable<Integer> {

    @Option(names = {"-i", "--input"}, description = "Input KCF file", required = true)
    private String inFile;

    @Option(names = {"-o", "--output"}, description = "Output KCF file", required = true)
    private String outFile;

    @Option(names = {"-w", "--window"}, description = "Window size", required = true)
    private int windowSize;

    @Option(names = {"-p", "--steps"}, description = "Step Size", defaultValue = "0")
    private int stepSize;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use", defaultValue = "1")
    private int nThreads;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    @Override
    public Integer call() throws Exception {
        try (KCFReader reader = new KCFReader(inFile)) {
            KCFHeader header = reader.getHeader();
            if (header.getStepSize() > 0){
                Logger.error(CLASS_NAME, "Cannot increase window size of a KCF file with overlapping windows (stepSize > 0)");
            }
            validateWindowSize(header.getWindowSize());
            Map<String, Window> newWindows = processWindows(reader, header);
            writeOutput(header, newWindows);
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error processing KCF file " + inFile);
            throw e;
        }
        return 0;
    }

    /***
     * Validate the window size, because we can only increase the window size not decrease it
     */
    private void validateWindowSize(int currentWindowSize) {
        if (currentWindowSize > windowSize) {
            Logger.error(CLASS_NAME, "Window size is smaller than the current window size");
        }
    }

    /***
     * Process the KCF file and increase the window size by merging windows
     */
    private Map<String, Window> processWindows(KCFReader reader, KCFHeader header) throws Exception {
        if (stepSize == 0){
            Map<String, Window> newWindows = new LinkedHashMap<>();
            LinkedList<Window> windowsBuffer = new LinkedList<>();

            int currentWindowSize = header.getWindowSize();
            int stepSize = windowSize / currentWindowSize + 1;
            int stepsReached = 0;
            int newWindowSize = 0;
            String lastSequenceName = null;

            for (Window window : reader) {
                String currentSequenceName = window.getSequenceName();

                if (lastSequenceName == null) {
                    lastSequenceName = currentSequenceName;
                }

                boolean sequenceChanged = !currentSequenceName.equals(lastSequenceName);
                if (sequenceChanged || stepsReached == stepSize) {
                    if (!windowsBuffer.isEmpty()) {
                        Window mergedWindow = combineWindows(windowsBuffer, header.getSamples(), header.getWeights());
                        newWindows.put(mergedWindow.getWindowId(), mergedWindow);
                        newWindowSize = Math.max(newWindowSize, mergedWindow.getEffLength());
                        windowsBuffer.clear();
                    }
                    stepsReached = 0;
                    lastSequenceName = currentSequenceName;
                }

                windowsBuffer.add(window);
                stepsReached++;
            }

            if (!windowsBuffer.isEmpty()) {
                Window mergedWindow = combineWindows(windowsBuffer, header.getSamples(), header.getWeights());
                newWindows.put(mergedWindow.getWindowId(), mergedWindow);
                newWindowSize = Math.max(newWindowSize, mergedWindow.getEffLength());
            }

            header.setWindowSize(newWindowSize);
            return newWindows;
        }
        else {

            ArrayList<Window> windowsBuffer = new ArrayList<>();

            int currentWindowSize = header.getWindowSize();
            String[] contigs = header.getContigs();
            // hashmap to hold windows per contig
            Map<String, ArrayList<Window>> contigWindows = new LinkedHashMap<>();
            for (Window window : reader) {
                contigWindows.computeIfAbsent(window.getSequenceName(), k -> new ArrayList<>()).add(window);
            }
            Logger.info(CLASS_NAME, "Merging windows");

            int windowsToMergeCount = windowSize / currentWindowSize + 1;
//        int stepCount = 0;
            int stepCount = Math.max(1, stepSize / currentWindowSize);
            int newWindowSize = 0;

            ExecutorService executor = Executors.newFixedThreadPool(nThreads);
            ConcurrentHashMap<String, ArrayList<Window>> newWindows = new ConcurrentHashMap<>();

            for (String contig : contigs) {
                int finalStepCount = stepCount;
                executor.submit(() -> {
                    ArrayList<Window> contigWindowList = contigWindows.get(contig);
                    if (contigWindowList == null || contigWindowList.isEmpty()) {
                        return; // no windows for this contig
                    }
                    for (int i = 0; i < contigWindowList.size(); i += finalStepCount) {
                        ArrayList<Window> windowsToMerge = new ArrayList<>();
                        for (int j = 0; j < windowsToMergeCount && (i + j) < contigWindowList.size(); j++) {
                            windowsToMerge.add(contigWindowList.get(i + j));
                        }
                        if (!windowsToMerge.isEmpty()) {
                            Window mergedWindow = combineWindows(windowsToMerge, header.getSamples(), header.getWeights());
                            newWindows.computeIfAbsent(contig, k -> new ArrayList<>()).add(mergedWindow);
                        }
                    }
                });
            }

            while (!executor.isTerminated()) {
                executor.shutdown();
                try {
                    if (!executor.awaitTermination(60, TimeUnit.SECONDS)) {
                        executor.shutdownNow();
                    }
                } catch (InterruptedException e) {
                    executor.shutdownNow();
                }
            }

            // free memory
            contigWindows.clear();
            windowsBuffer.clear();

            Logger.info(CLASS_NAME, "Compiling final windows");
            int newStepSize = 0;
            int lastStart = 0;
            int lastEffLen = 0;
            Map<String, Window> finalNewWindows = new LinkedHashMap<>();
            for (String contig : contigs) {
                ArrayList<Window> contigNewWindows = newWindows.get(contig);
                if (contigNewWindows != null) {
                    for (Window window : contigNewWindows) {
                        finalNewWindows.put(window.getWindowId(), window);
                        if (window.getEffLength() > newWindowSize) {
                            newWindowSize = window.getEffLength();
                        }
                        if (window.getStart() > lastStart || window.getStart() > 0) {
                            int tmpSize = window.getStart() - lastStart;
                            if (tmpSize > newStepSize) {
                                newStepSize = tmpSize;
                            }
                        }
                        lastStart = window.getStart();
                        lastEffLen = window.getEffLength();
                    }
                }
            }
            header.setWindowSize(newWindowSize);
            header.setStepSize(newStepSize);
            Logger.info(CLASS_NAME, "New window size: " + newWindowSize + ", New step size: " + newStepSize);
            return finalNewWindows;
        }
    }

    /***
     * Write the output KCF file with the new windows
     */
    private void writeOutput(KCFHeader header, Map<String, Window> windows) throws Exception {
        try (KCFWriter writer = new KCFWriter(outFile)) {
            header.addCommandLine(HelperFunctions.getCommandLine());
            writer.writeHeader(header);

            String[] headerSamples = header.getSamples();
            for (Window window : windows.values()) {
                window.alignSamplesWithHeader(headerSamples);
            }

            writer.writeWindows(windows.values());
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Error writing KCF file: " + outFile);
            throw e;
        }
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /***
     * Combine multiple windows into a single window
     */
    private Window combineWindows(
            List<Window> windows,
            String[] headerSamples,
            double[] weights) {

        if (windows.isEmpty()) {
            throw new IllegalArgumentException("Window list cannot be empty");
        }

        final int sampleCount = headerSamples.length;
        final int totalWindows = windows.size();

        final Window first = windows.get(0);
        final Window last = windows.get(totalWindows - 1);

        int newStart = first.getStart();
        int newEnd   = last.getEnd();

        int tot = 0;
        int[] va = new int[sampleCount];
        int[] ob = new int[sampleCount];
        int[] id = new int[sampleCount];
        int[] ld = new int[sampleCount];
        int[] rd = new int[sampleCount];
        long[] kt = new long[sampleCount];
        int[] prevRightDistance = new int[sampleCount];

        boolean singleWindow = (totalWindows == 1);

        for (int wIndex = 0; wIndex < totalWindows; wIndex++) {
            Window w = windows.get(wIndex);
            tot += w.getTotalKmers();

            boolean isFirst = (wIndex == 0);
            boolean isLast  = (wIndex == totalWindows - 1);

            for (int i = 0; i < sampleCount; i++) {
                Data data = w.getData().get(headerSamples[i]);  // ideally O(1) array lookup
                if (data == null) continue;

                int left  = data.getLeftDistance();
                int right = data.getRightDistance();
                int vars  = data.getVariations();

                if (prevRightDistance[i] > 0 && left > 0 && vars > 0) {
                    va[i] += vars - 1;
                } else {
                    va[i] += vars;
                }

                ob[i] += data.getObservedKmers();
                id[i] += data.getInnerDistance();
                kt[i] += (long) data.getMeanKmerCount() * data.getObservedKmers();

                if (singleWindow) {
                    ld[i] += left;
                    rd[i] += right;
                } else if (isFirst) {
                    ld[i] += left;
                    id[i] += right;
                } else if (isLast) {
                    rd[i] += right;
                    id[i] += left;
                } else {
                    id[i] += left + right;
                }

                prevRightDistance[i] = right;
            }
        }

        Window newWindow = new Window(
                first.getSequenceName() + "_" + newStart,
                first.getSequenceName(),
                newStart,
                newEnd
        );
        newWindow.addTotalKmers(tot);
        newWindow.setEffLength(newEnd - newStart);

        for (int i = 0; i < sampleCount; i++) {
            newWindow.addData(
                    headerSamples[i],
                    ob[i],
                    va[i],
                    id[i],
                    ld[i],
                    rd[i],
                    kt[i],
                    "N",
                    weights
            );
        }
        return newWindow;
    }
}
// EOF