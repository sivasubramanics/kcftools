package nl.wur.bis.kcftools.Plugins;


import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.*;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

@Command(name = "kcf2plink", description = "Convert KCF windows to PED format")
public class KCFToPed implements Runnable, Callable<Integer> {
    @Option(names = {"-i", "--input"}, description = "Input KCF file", required = true)
    private String inFile;
    @Option(names = {"-o", "--output"}, description = "Output PED file prefix", required = true)
    private String outPrefix;
    @Option(names = {"-a", "--score_a"}, description = "Lower score cut-off for reference allele", required = false)
    private double scoreA = 95.0;
    @Option(names = {"-b", "--score_b"}, description = "Lower score cut-off for alternate allele", required = false)
    private double scoreB = 60.0;
    @Option(names = {"--score_n"}, description = "Score value for missing data (default = 30.0)", required = false)
    private double scoreN = 30.0;
    @Option(names = {"--chrs"}, description = "List file with chromosomes to include", required = false)
    private String chrsFile = null;
    @Option(names = {"--maf"}, description = "Minimum allele frequency to consider a window valid", required = false)
    private double minMAF = 0.05;
    @Option(names = {"--max-missing"}, description = "Maximum proportion of missing data to consider a window valid", required = false)
    private double maxMissing = 0.8;

    private final String CLASSNAME = this.getClass().getSimpleName();
    public KCFToPed() {
    }

    public KCFToPed(String inFile, String outPrefix) {
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
        convertKCFtoPED(inFile, outPrefix, scoreA, scoreB);
        return 0;
    }

    private void convertKCFtoPED(String inFile, String outPrefix, double scoreA, double scoreB) {
        Logger.warning(CLASSNAME, "This is an experimental feature, use with caution!");
        try (KCFReader reader = new KCFReader(inFile);
             BufferedWriter writer = new BufferedWriter(new FileWriter(outPrefix + ".ped"));
             BufferedWriter mapWriter = new BufferedWriter(new FileWriter(outPrefix + ".map"));
                BufferedWriter contigsMapWriter = new BufferedWriter(new FileWriter(outPrefix + ".contigsMap"))) {

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
                map[i] = contigID + "\t" + i + "\t" + "0" + "\t" + window.getStart();
                if (!contigsMap.contains(window.getSequenceName() + "\t" + contigID)) {
                    contigsMap.add(window.getSequenceName() + "\t" + contigID);
                }
                int[] alleles = new int[sampleCount];
                for (int j = 0; j < sampleCount; j++) {
                    double score = window.getData().get(samples[j]).getScore();
                    if (score >= scoreA) {
                        matrix[j][i] = 0;
                        alleles[j] = 0; // reference allele 100-scoreA
                    } else if (score >= scoreB) {
                        matrix[j][i] = 2;
                        alleles[j] = 2; // alternate allele scoreB-scoreA
                    } else if (score <= scoreN) {
                        matrix[j][i] = -1;
                        alleles[j] = -1; // missing data 0-scoreN
                    } else {
                        matrix[j][i] = 1;
                        alleles[j] = 1; // heterozygous scoreN-scoreB
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
            for (int j = 0; j < sampleCount; j++) {
                writer.write(samples[j] + "\t" + samples[j] + "\t0\t0\t0\t-9"); // PED format requires two identical columns for the sample name
                for (int k = 0; k < i; k++) {
                    if (!badWindows.contains(k)) {
                        switch (matrix[j][k]) {
                        case 0:
                            writer.write("\tA\tA"); // reference allele
                            break;
                        case 2:
                            writer.write("\tG\tG"); // alternate allele
                            break;
                            case 1:
                            writer.write("\tA\tG"); // missing data
                            break;
                        case -1:
                            writer.write("\t0\t0"); // missing data
                            break;
                        default:
                            writer.write("\t0\t0"); // default case for unexpected values
                            break;
                        }
                    }
                }
                writer.newLine();
            }
            writer.flush();
            Logger.info(CLASSNAME, "Generated Matrix file: " + outPrefix + ".matrix.tsv");

        } catch (Exception e) {
            e.printStackTrace();
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
}
