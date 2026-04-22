package nl.wur.bis.kcftools.Plugins;


import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.KCFWriter;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.util.*;
import java.util.concurrent.Callable;

@Command(name = "extractKCF", description = "Extract KCF for one or list of samples")
public class ExtractKCF implements Callable<Integer>, Runnable {
    // input KCF file
    @Option(names = {"-k", "--kcf"}, description = "Input KCF file", required = true)
    private String kcfFile;
    // comma-separated list of sample names to extract
    @Option(names = {"-s", "--samples"}, description = "Comma-separated list of sample names to extract", required = false)
    private String sampleNames;
    // input list file with sample names and corresponding fasta files
    @Option(names = {"-l", "--list"}, description = "Input list file with sample names and corresponding fasta files", required = false)
    private String sampleListFile;
    // output KCF file
    @Option(names = {"-o", "--output"}, description = "Output KCF file", required = true)
    private String outputKCF;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public ExtractKCF() {
    }

    public ExtractKCF(String kcfFile, String sampleNames, String sampleListFile, String outputKCF) {
        this.kcfFile = kcfFile;
        this.sampleNames = sampleNames;
        this.sampleListFile = sampleListFile;
        this.outputKCF = outputKCF;
        run();
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public Integer call() throws Exception {
        extractKCF();
        return 0;
    }

    private void extractKCF() {
        // check if sampleNames or sampleListFile is provided, and make sure only one is provided
        if ((sampleNames == null || sampleNames.isEmpty()) && (sampleListFile == null || sampleListFile.isEmpty())) {
            Logger.error(CLASS_NAME, "Either sample names or sample list file must be provided.");
        }

        // if sampleNames is provided, split by comma and trim spaces
        Set<String> samplesToExtract = new HashSet<>();
        if (sampleNames != null && !sampleNames.isEmpty()) {
            String[] samplesArray = sampleNames.split(",");
            for (String sample : samplesArray) {
                samplesToExtract.add(sample.trim());
            }
        }

        // if sampleListFile is provided, read the file and extract sample names
        if (sampleListFile != null && !sampleListFile.isEmpty()) {
            try (java.io.BufferedReader br = new java.io.BufferedReader(new java.io.FileReader(sampleListFile))) {
                java.util.List<String> sampleList = new java.util.ArrayList<>();
                String line;
                while ((line = br.readLine()) != null) {
                    String[] parts = line.split("\t");
                    if (parts.length > 0) {
                        sampleList.add(parts[0].trim());
                    }
                }
                // if it's already in samplesToExtract, log warning and skip
                for (String sample : sampleList) {
                    if (samplesToExtract.contains(sample)) {
                        Logger.warning(CLASS_NAME, "Sample " + sample + " is already in the samples to extract list, skipping.");
                    } else {
                        samplesToExtract.add(sample);
                    }
                }
            } catch (java.io.IOException e) {
                Logger.error(CLASS_NAME, "Error reading sample list file: " + e.getMessage());
            }
        }

        // read the KCF file and extract the relevant samples
        try (KCFReader kcfReader = new KCFReader(kcfFile)) {
            KCFHeader header = kcfReader.getHeader();
            // read through the samples array from header, and if the sample is in samplesToExtract, consider to extract, else throw warning
            List<String> headerSamples = Arrays.asList(header.getSamples());
            Iterator<String> iterator = samplesToExtract.iterator();
            while (iterator.hasNext()) {
                String sample = iterator.next();
                if (!headerSamples.contains(sample)) {
                    Logger.warning(CLASS_NAME, "Sample " + sample + " not found in KCF file header, skipping.");
                    iterator.remove();
                }
            }
            if (samplesToExtract.isEmpty()) {
                Logger.error(CLASS_NAME, "No valid samples to extract found in KCF file.");
                return;
            }
            KCFHeader newHeader = new KCFHeader(header);
            newHeader.setSamples(samplesToExtract.toArray(new String[0]));
            Logger.info(CLASS_NAME, "Extracting " + samplesToExtract.size() + " samples to " + outputKCF);
            try (KCFWriter kcfWriter = new KCFWriter(outputKCF)) {
                newHeader.addCommandLine(HelperFunctions.getCommandLine());
                kcfWriter.writeHeader(newHeader);
                for (Window window : kcfReader) {
                    window.alignSamplesWithHeader(newHeader.getSamples());
                    kcfWriter.writeWindow(window);
                }
            } catch (Exception e) {
                Logger.error(CLASS_NAME, "Error writing output KCF file: " + e.getMessage());
            }
        }
        catch (Exception e) {
            Logger.error(CLASS_NAME, "Error reading KCF file: " + e.getMessage());
        }
    }
}
// EOF