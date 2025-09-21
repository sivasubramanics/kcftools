package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

/***
 * Convert the KCF windows to genotype matrix
 */
@Command(name = "kcfToMatrix", description = "Convert the KCF windows to genotype matrix")
public class KCFToMatrix implements Callable<Integer>, Runnable {

    @Option(names = {"-i", "--input"}, description = "Input KCF file", required = true)
    private String inFile;

    @Option(names = {"-o", "--output"}, description = "Output file prefix", required = true)
    private String outPrefix;

    @Option(names = {"-a", "--score_a"}, description = "Lower score cut-off for reference allele (default = 95.0)", required = false)
    private double scoreA = 95.0;

    @Option(names = {"-b", "--score_b"}, description = "Lower score cut-off for alternate allele (default = 60.0)", required = false)
    private double scoreB = 60.0;

    @Option(names = {"--score_n"}, description = "Score value for missing data (default = 30.0)", required = false)
    private double scoreN = 30.0;

    @Option(names = {"-r", "--rdata"}, description = "Convert the matrix to RData", required = false)
    private boolean rdata = false;

    @Option(names = {"--maf"}, description = "minimum allele frequency to consider a window valid", required = false)
    private double minMAF = 0.05;

    @Option(names = {"--max-missing"}, description = "maximum proportion of missing data to consider a window valid", required = false)
    private double maxMissing = 0.8;

    @Option(names = {"--chrs"}, description = "List file with chromosomes to include", required = false)
    private String chrsFile = null;

    private final String CLASSNAME = this.getClass().getSimpleName();

    public KCFToMatrix() {
    }

    public KCFToMatrix(String inFile, String outPrefix) {
        this.inFile = inFile;
        this.outPrefix = outPrefix;
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
        convertKCFToMatrix();
        return 0;
    }

    /***
     * Main method to convert KCF to matrix
     */
    private void convertKCFToMatrix() {
        try (KCFReader reader = new KCFReader(inFile);
             BufferedWriter writer = new BufferedWriter(new FileWriter(outPrefix + ".matrix.tsv"));
             BufferedWriter mapWriter = new BufferedWriter(new FileWriter(outPrefix + ".map.tsv"));
             BufferedWriter contigsMapWriter = new BufferedWriter(new FileWriter(outPrefix + ".contigsMap.tsv"))) {

            KCFHeader header = reader.getHeader();
            int sampleCount = header.getSamples().length;
            int windowCount = header.getWindowCount();
            int[][] matrix = new int[sampleCount][windowCount];
            String[] map = new String[windowCount];
            List<String> contigsMap = new LinkedList<>();
            Set<Integer> badWindows = new HashSet<>();
            String[] samples = header.getSamples();
            int contigID;
            int i = 0;
            Set<String> chrs = null;

            if (chrsFile != null) {
                try (BufferedReader br = new BufferedReader(new FileReader(chrsFile))) {
                    chrs = new HashSet<>();
                    String line;
                    while ((line = br.readLine()) != null) {
                        if (line.startsWith("#") || line.trim().isEmpty()) {
                            continue; // skip comments and empty lines
                        }
                        chrs.add(line.trim());
                    }
                }
            }

            for (Window window : reader) {
                contigID = header.getContigID(window.getSequenceName()) + 1;
                map[i] = i + "\t" + contigID + "\t" + window.getStart();
                if (!contigsMap.contains(window.getSequenceName() + "\t" + contigID)) {
                    contigsMap.add(window.getSequenceName() + "\t" + contigID);
                }
                int[] alleles = new int[sampleCount];
                for (int j = 0; j < sampleCount; j++) {
                    double score = window.getData().get(samples[j]).getScore();
                    if (score >= scoreA) {
                        matrix[j][i] = 0;
                        alleles[j] = 0;
                    } else if (score >= scoreB) {
                        matrix[j][i] = 2;
                        alleles[j] = 2;
                    } else if (score <= scoreN) {
                        matrix[j][i] = -1;
                        alleles[j] = -1; // missing data
                    } else {
                        matrix[j][i] = 1;
                        alleles[j] = 1;
                    }
                }
                if (chrs != null && !chrs.contains(window.getSequenceName())) {
                    badWindows.add(i);
                    continue; // skip this window if not in the specified chromosomes
                }
                if (badWindow(alleles)) {
                    badWindows.add(i);
                }
                i++;
            }

            // Write the map file
            mapWriter.write("name\tchromosome\tposition");
            mapWriter.newLine();
            for (int m = 0; m < i; m++) {
                if (map[m] != null && !badWindows.contains(m)) {
                    mapWriter.write(map[m]);
                    mapWriter.newLine();
                }
            }
            mapWriter.flush();
            Logger.info(CLASSNAME, "Generated Map file: " + outPrefix + ".map.tsv");


            // Write the contigs map file
            for (String contig : contigsMap) {
                contigsMapWriter.write(contig);
                contigsMapWriter.newLine();
            }
            contigsMapWriter.flush();
            Logger.info(CLASSNAME, "Generated Contigs Map file: " + outPrefix + ".contigsMap.tsv");

            // Write the matrix file
            writer.write("taxa");
            for (int k = 0; k < i; k++) {
                if (!badWindows.contains(k)) {
                    writer.write("\t" + k);
                }
            }
            writer.newLine();
            for (int j = 0; j < sampleCount; j++) {
                writer.write(samples[j]);
                for (int k = 0; k < i; k++) {
                    if (!badWindows.contains(k)) {
                        writer.write("\t" + (matrix[j][k] == -1 ? 1 : matrix[j][k]));
                    }
                }
                writer.newLine();
            }
            writer.flush();
            Logger.info(CLASSNAME, "Generated Matrix file: " + outPrefix + ".matrix.tsv");

        } catch (Exception e) {
            e.printStackTrace();
        }

        if (rdata) {
            convertGTmatrixToRdata(outPrefix + ".matrix.tsv", outPrefix + ".map.tsv");
        }
    }

    private boolean badWindow(int[] alleles) {
        int count0 = 0, count1 = 0, count2 = 0, countN = 0;
        for (int allele : alleles) {
            if (allele == 0) count0++;
            else if (allele == 1) count1++;
            else if (allele == 2) count2++;
            else if (allele == -1) countN++;
        }

        int validCount = alleles.length - countN;
        return (count0 == alleles.length || count1 == alleles.length || count2 == alleles.length || countN == alleles.length)
                || (validCount > 0 && (count0 <= minMAF * validCount || count2 <= minMAF * validCount))
                || (countN >= maxMissing * alleles.length || (countN + count1) >= maxMissing * alleles.length);
    }

    private void convertGTmatrixToRdata(String matrixFile, String mapFile) {
        if (!HelperFunctions.isInstalled("Rscript")) {
            Logger.error(CLASSNAME, "Rscript is not installed. Please install Rscript and try again.");
            return;
        }
        Logger.info(CLASSNAME, "Converting matrix to RData");

        String rscriptName = "convertGTmatrixToRdata_" + System.currentTimeMillis() + ".R";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(rscriptName))) {
            writer.write("df <- read.table(\"" + matrixFile + "\", head = TRUE, sep = \"\\t\")\n");
            writer.write("save(df, file = \"" + matrixFile.replaceAll("\\.tsv$", ".RData") + "\")\n");
            writer.write("df <- read.table(\"" + mapFile + "\", head = TRUE, sep = \"\\t\")\n");
            writer.write("save(df, file = \"" + mapFile.replaceAll("\\.tsv$", ".RData") + "\")\n");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        try {
            Logger.info(CLASSNAME, "Running Rscript: " + rscriptName);
            HelperFunctions.tryExec("Rscript " + rscriptName);
            HelperFunctions.deleteFile(rscriptName);
        } catch (Exception e) {
            Logger.error(CLASSNAME, "Error while running Rscript: " + e.getMessage());
        }
    }
}
// EOF
