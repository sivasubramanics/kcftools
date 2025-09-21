package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;

/***
 * Find IBS windows in a KCF file
 */
@Command(name = "findIBS", description = "Find IBS windows in a KCF file")
public class FindIBS implements Callable<Integer>, Runnable {
    // in KCF file name
    @Option(names = {"-i", "--input"}, description = "Input KCF file name", required = true)
    private String inFile;
    // in output file name
    @Option(names = {"-o", "--output"}, description = "Output KCF file name", required = true)
    private String outFile;
    // flag to detect IBS or Variable Regions
    @Option(names = {"--var"}, description = "Detect Variable Regions instead of IBS [default: false]", required = false, defaultValue = "false")
    private boolean detectVar;
    // number of minimum consecutive windows
    @Option(names = {"--min"}, description = "Minimum number of consecutive windows [default: 4]", required = false, defaultValue = "4")
    private int minConsecutive;
    // score cut-off
    @Option(names = {"--score"}, description = "Score cut-off [default: 95.00]", required = false, defaultValue = "95")
    private float scoreCutOff;
    // flag to write summary tsv file
    @Option(names = {"--summary"}, description = "Write summary tsv file [default: false]", required = false, defaultValue = "false")
    private boolean writeSummary;
    // flag to write bed file
    @Option(names = {"--bed"}, description = "Write bed file [default: false]", required = false, defaultValue = "false")
    private boolean writeBed;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public FindIBS() {
    }

    public FindIBS(String inFile, String outFile, boolean detectVar, int minConsecutive, float scoreCutOff) {
        this.inFile = inFile;
        this.outFile = outFile;
        this.detectVar = detectVar;
        this.minConsecutive = minConsecutive;
        this.scoreCutOff = scoreCutOff;
    }

    @Override
    public Integer call() throws Exception {
        findIBS();
        return 0;
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
     * Main function to find IBS windows in a KCF file
     */
    private void findIBS() throws Exception {

        // if outFile doesnt end with .kcf, add it
        if (!outFile.endsWith(".kcf")) {
            outFile += ".kcf";
        }

        KCFReader reader = new KCFReader(inFile);
        KCFHeader header = reader.getHeader();
        if (header.getStepSize() > 0){
            minConsecutive = header.getWindowSize() / header.getStepSize();
            Logger.warning(CLASS_NAME, "Input KCF file is created with step size. Hence we are using the --min = windowSize/stepSize [" + minConsecutive + "]");
        }
        HashMap<String, Integer> windowsCount = new HashMap<>();
        HashMap<String, ArrayList<Window>> windowsList = new HashMap<>();
        for (Window window: reader) {
            String name = window.getSequenceName();
            if (!windowsCount.containsKey(name)) {
                windowsCount.put(name, 0);
                windowsList.put(name, new ArrayList<>());
            }
            windowsCount.put(name, windowsCount.get(name) + 1);
            windowsList.get(name).add(window);
        }
        reader.close();

        // create a Array of Window objects for each chromosome
        HashMap<String, Window[]> windows = new HashMap<>();
        for (String name : windowsCount.keySet()) {
            int windowSize = windowsCount.get(name);
            windows.put(name, new Window[windowSize]);
        }


        for (String chromName : windowsList.keySet()) {
            int winNum = 0;
            for (Window window: windowsList.get(chromName)) {
                if (window == null) {
                    Logger.error(CLASS_NAME, "Window not found in index: " + chromName);
                    continue;
                }
                windows.get(chromName)[winNum] = window;
                winNum++;
            }
        }

        String[] samples = header.getSamples();
        for (String sample : samples) {
            Logger.info(CLASS_NAME, "Finding IBS for sample: " + sample);
            int blockNum = 0;
            String blockChrom = null;
            boolean firstIBSFound = false;
            for (String chromName : windows.keySet()) {
                Window[] chromWindows = windows.get(chromName);
                int numNA = 0;
                double score;
                boolean isIBSRegion;

                for (int i = 0; i < chromWindows.length; i++) {
                    Window window = chromWindows[i];

                    if (window == null) {
                        Logger.error(CLASS_NAME, "Window not found in index: " + chromName + ":" + i);
                        continue;
                    }

                    if (!window.getData().containsKey(sample)) {
                        continue;
                    }

                    score = window.getScore(sample);
                    isIBSRegion = detectVar ? score < scoreCutOff  : score >= scoreCutOff;

                    if (isIBSRegion) {
                        if (!firstIBSFound) {
                            blockNum = 1;
                            firstIBSFound = true;
                        } else if (numNA >= minConsecutive || (blockChrom != null && !blockChrom.equals(chromName))) {
                            blockNum++;
                        }
                        blockChrom = chromName;
                        window.setIBS(sample, blockNum);
                        numNA = 0;
                    } else {
                        numNA++;
                        window.setIBS(sample, -1);
                    }
                }
            }
        }


        try (KCFWriter writer = new KCFWriter(outFile)) {
            header.setIBS(true);
            header.addCommandLine(HelperFunctions.getCommandLine());
            writer.writeHeader(header);
            for (String chromName : windows.keySet()) {
                for (Window window : windows.get(chromName)) {
                    writer.writeWindow(window);
                }
            }
        }

        if (writeSummary) {
            try (BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(outFile.replace(".kcf", ".summary.tsv")))) {
                // write header for the summary file
                summaryWriter.write("Block\tSample\tChromosome\tStart\tEnd\tLength\tTotalBlocks\tIBSBlocks\tIBSProportion\tMeanScore\n");

                for (String sample : samples) {
                    Map<Integer, List<Window>> blocks = new LinkedHashMap<>();
                    for (String chromName : windows.keySet()) {
                        Window[] chromWindows = windows.get(chromName);
                        List<Window> naWindows = new ArrayList<>();

                        for (Window window : chromWindows) {
                            int ibsValue = window.getIbs(sample);

                            if (ibsValue == -1) {
                                naWindows.add(window);
                            } else {
                                if (blocks.containsKey(ibsValue)) {
                                    blocks.get(ibsValue).addAll(naWindows);
                                    blocks.get(ibsValue).add(window);
                                    naWindows.clear();
                                } else {
                                    blocks.put(ibsValue, new ArrayList<>());
                                    blocks.get(ibsValue).add(window);
                                    naWindows.clear();
                                }
                            }
                        }
                    }

                    // write BED file if required
                    if (writeBed) {
                        writeBedFile(outFile, sample, blocks);
                    }

                    // write summary for each block
                    for (Map.Entry<Integer, List<Window>> entry : blocks.entrySet()) {
                        writeSummaryEntry(summaryWriter, entry, sample);
                    }
                }
            }
        }
    }

    /***
     * Write the BED file for the IBS blocks
     */
    private void writeBedFile(String outFile, String sample, Map<Integer, List<Window>> blocks) throws IOException {
        try (BufferedWriter bedWriter = new BufferedWriter(new FileWriter(outFile.replace(".kcf", "." + sample + ".bed")))) {
            for (List<Window> block : blocks.values()) {
                if (!block.isEmpty()) {
                    int start = block.get(0).getStart();
                    int end = block.get(block.size() - 1).getEnd();
                    bedWriter.write(block.get(0).getSequenceName() + "\t" + start + "\t" + end + "\n");
                }
            }
        }
    }

    /***
     * Write the summary entry for each block
     */
    private void writeSummaryEntry(BufferedWriter summaryWriter, Map.Entry<Integer, List<Window>> entry, String sample) throws IOException {
        List<Window> block = entry.getValue();
        if (block.isEmpty()) return;

        int totalBlocks = block.size();
        int ibsBlocks = 0;
        float meanScore = 0;

        int start = block.get(0).getStart();
        int end = block.get(block.size() - 1).getEnd();

        for (Window window : block) {
            meanScore += window.getScore(sample);
            if (window.getIbs(sample) != -1) {
                ibsBlocks++;
            }
        }

        meanScore /= totalBlocks; // Calculate mean score
        float ibsProportion = (float) ibsBlocks / totalBlocks;

        // Write summary entry
        summaryWriter.write(String.format(
                "%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n",
                entry.getKey(),
                sample,
                block.get(0).getSequenceName(),
                start,
                end,
                (end - start),
                totalBlocks,
                ibsBlocks,
                ibsProportion,
                meanScore
        ));
    }
}
//EOF