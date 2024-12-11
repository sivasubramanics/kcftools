package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.KCFWriter;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import picocli.CommandLine.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.Callable;

@Command(name = "cohort", description = "Create a cohort of samples kcf files")
public class Cohort implements Callable<Integer>, Runnable {
    @Option(names = {"-o", "--output"}, description = "Output file name", required = true)
    private String outFile;

    @Option(names = {"-i", "--input"}, description = "List of samples kcf files", split = ",", required = false)
    private String[] inFiles;

    @Option(names = {"-l", "--list"}, description = "File containing list of samples kcf files", required = false)
    private String listFile;

    private static final String CLASS_NAME = Cohort.class.getSimpleName();

    @Override
    public Integer call() throws Exception {
        if (inFiles == null && listFile == null) {
            HelperFunctions.log("error", CLASS_NAME, "No input files provided");
        }

        if (listFile != null) {
            inFiles = readListFile(listFile);
        }

        cohortKcfFiles();
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

    private String[] readListFile(String listFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(listFile))) {
            return reader.lines().toArray(String[]::new);
        }
    }

    private void cohortKcfFiles() {
        Map<String, Window> windows = new LinkedHashMap<>();
        KCFHeader header = null;
        KCFHeader tmpHeader = null;

        for (int i = 0; i < inFiles.length; i++) {
            try (KCFReader reader = new KCFReader(inFiles[i])) {
                if (i == 0) {
                    header = reader.getHeader();
                    for (Window window : reader) {
                        windows.put(getWindowKey(window), window);
                    }
                } else {
                    tmpHeader = reader.getHeader();
                    assert header != null;
                    if(!header.equals(tmpHeader)) {
                        HelperFunctions.log("error", CLASS_NAME, "Headers mismatch found in sample: " + inFiles[i]);
                    }
                    header.mergeHeader(tmpHeader);
                    for (Window window : reader) {
                        String key = getWindowKey(window);
                        if (windows.containsKey(key)) {
                            windows.get(key).addData(window.getData());
                        } else {
                            HelperFunctions.log("error", CLASS_NAME, "Windows mismatch found in sample: " + inFiles[i] + " at window: " + window);
                        }
                    }
                }
            } catch (Exception e) {
                HelperFunctions.log("error", CLASS_NAME, "Error reading KCF file: " + inFiles[i]);
            }
        }

//        // reshuffle samples in windows to match the order of samples in header
//        for (Window window : windows.values()) {
//            window.shuffleSamples(header.getSamples());
//        }

        try(KCFWriter writer = new KCFWriter(outFile)) {
            assert header != null;
            header.addCommandLine(HelperFunctions.getCommandLine());
            writer.writeHeader(header);

            // Align samples in each window with the header
            String[] headerSamples = header.getSamples();
            for (Window window : windows.values()) {
                window.alignSamplesWithHeader(headerSamples);
            }

            writer.writeWindows(windows.values());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private String getWindowKey(Window window) {
        return window.getSequenceName() + ":" + window.getStart() + "-" + window.getEnd();
    }

}