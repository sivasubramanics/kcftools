package nl.wur.bis.kcftools.Plugins;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.concurrent.Callable;
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

    private final String CLASS_NAME = this.getClass().getSimpleName();

    @Override
    public Integer call() throws Exception {
        try (KCFReader reader = new KCFReader(inFile)) {
            KCFHeader header = reader.getHeader();
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
    private Window combineWindows(LinkedList<Window> windows, String[] headerSamples, double[] weights) {
        if (windows.isEmpty()) {
            throw new IllegalArgumentException("Window list cannot be empty");
        }

        int newStart = windows.getFirst().getStart();
        int newEnd = windows.getLast().getEnd();
        int totalWindows = windows.size();

        int tot = 0;
        int[] va = new int[headerSamples.length];
        int[] ob = new int[headerSamples.length];
        int[] id = new int[headerSamples.length];
        int[] ld = new int[headerSamples.length];
        int[] rd = new int[headerSamples.length];

        boolean singleWindow = (totalWindows == 1);
        int winIndex = 0;
        int[] prevRightDistance = new int[headerSamples.length]; // right distance of previous window

        for (Window window : windows) {
            tot += window.getTotalKmers();

            for (int i = 0; i < headerSamples.length; i++) {
                Data data = window.getData().get(headerSamples[i]);
                if (data == null) {
                    continue;
                }

                int left = data.getLeftDistance();
                int right = data.getRightDistance();
                int vars = data.getVariations();

                // adjust variation count for overlap
                if (prevRightDistance[i] > 0 && left > 0 && vars > 0) {
                    va[i] += (vars - 1); // subtract 1 for overlap
                } else {
                    va[i] += vars;
                }

                ob[i] += data.getObservedKmers();
                id[i] += data.getInnerDistance();

                if (singleWindow) {
                    ld[i] += left;
                    rd[i] += right;
                } else if (winIndex == 0) {
                    ld[i] += left;
                    id[i] += right;
                } else if (winIndex == totalWindows - 1) {
                    rd[i] += right;
                    id[i] += left;
                } else {
                    id[i] += left + right;
                }

                prevRightDistance[i] = right; // track for next iteration
            }
            winIndex++;
        }
        Window newWindow = new Window(windows.getFirst().getSequenceName() + "_" + newStart, windows.getFirst().getSequenceName(), newStart, newEnd);
        newWindow.addTotalKmers(tot);
        newWindow.setEffLength(newEnd - newStart);
        for (int i = 0; i < headerSamples.length; i++) {
            newWindow.addData(
                headerSamples[i],
                ob[i],
                va[i],
                id[i],
                ld[i],
                rd[i],
                "N",
                weights
            );
        }
        return newWindow;
    }
}
// EOF