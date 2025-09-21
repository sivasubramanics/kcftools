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

@Command(name = "kcfToGenotypeTable", description = "Convert KCF to Genotype Table")
public class KCFToGenotypeTable implements Callable<Integer>, Runnable {

    @Option(names = {"-i", "--input"}, description = "Input KCF file", required = true)
    private String inFile;

    @Option(names = {"-o", "--output"}, description = "Output file", required = true)
    private String outFile;

    @Option(names = {"--score_a"}, description = "Lower score cut-off for reference allele (default = 95.0)", required = false)
    private double scoreA = 95.0;

    @Option(names = {"--score_b"}, description = "Lower score cut-off for alternate allele (default = 60.0)", required = false)
    private double scoreB = 60.0;

    @Option(names = {"--score_n"}, description = "Score value for missing data (default = 30.0)", required = false)
    private double scoreN = 30.0;

    @Option(names = {"--maf"}, description = "minimum allele frequency to consider a window valid", required = false)
    private double minMAF = 0.00;

    @Option(names = {"--max-missing"}, description = "maximum proportion of missing data to consider a window valid", required = false)
    private double maxMissing = 1.00;

    @Option(names = {"--chrs"}, description = "List file with chromosomes to include", required = false)
    private String chrsFile = null;

    private final String CLASSNAME = this.getClass().getSimpleName();

    public KCFToGenotypeTable() {
    }

    public KCFToGenotypeTable(String inFile, String outFile) {
        this.inFile = inFile;
        this.outFile = outFile;
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
        convertKCFtoGenotypeTable();
        return 0;
    }

    private void convertKCFtoGenotypeTable() {

        validateScores();

        try (KCFReader reader = new KCFReader(inFile);
             BufferedWriter writer = new BufferedWriter(new java.io.FileWriter(outFile));
             BufferedWriter contigsMapWriter = new BufferedWriter(new FileWriter(outFile + ".contigsMap.tsv"))) {
            KCFHeader header = reader.getHeader();
            int sampleCount = header.getSamples().length;
            int windowCount = header.getWindowCount();
            String[] samples = header.getSamples();
            List<String> contigsMap = new LinkedList<>();
            int contigID;
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
            writer.write("# Genotype Table " +
                    "0:" + scoreA + " - 100.00, " +
                    "2:" + scoreB + " - " + scoreA + ", " +
                    "1:" + scoreN + " - " + scoreB + ", " +
                    "-1: <=" + scoreN);
            writer.newLine();
            writer.write("ID\tCHR\tSTART\tEND");
            for (String sample : samples) {
                writer.write("\t" + sample);
            }
            writer.newLine();

            for (Window window : reader) {
                contigID = header.getContigID(window.getSequenceName()) + 1;
                if (!contigsMap.contains(window.getSequenceName() + "\t" + contigID)) {
                    contigsMap.add(window.getSequenceName() + "\t" + contigID);
                }
                StringBuilder outLine = new StringBuilder();
                if (chrs != null && !chrs.contains(window.getSequenceName())) {
                    continue;
                }
                int[] alleles = new int[sampleCount];
                for (int j = 0; j < sampleCount; j++) {
                    double score = window.getData().get(samples[j]).getScore();
                    if (score >= scoreA) {
                        alleles[j] = 0; // homozygous ref (0) 100 - scoreA
                    } else if (score >= scoreB) {
                        alleles[j] = 2; // homozygous alt (2) scoreB - scoreA
                    } else if (score <= scoreN) {
                        alleles[j] = -1; // missing data (-1) 0 - scoreN
                    } else {
                        alleles[j] = 1; // heterozygous (1) scoreN - scoreB
                    }
                }
                if (badWindow(alleles) && (minMAF > 0.0 || maxMissing < 1.0)) {
                    continue;
                }
                outLine.append(window.getWindowId());
                outLine.append("\t").append(contigID);
                outLine.append("\t").append(window.getStart());
                outLine.append("\t").append(window.getEnd());
                for (int allele : alleles) {
                    outLine.append("\t").append(allele);
                }
                writer.write(outLine.toString());
                writer.newLine();
            }
            writer.flush();
            Logger.info(CLASSNAME, "Genotype table written to: " + outFile);

            contigsMapWriter.write("contigName\tcontigID");
            contigsMapWriter.newLine();
            for (String contig : contigsMap) {
                contigsMapWriter.write(contig);
                contigsMapWriter.newLine();
            }
            contigsMapWriter.flush();
            Logger.info(CLASSNAME, "Generated Contigs Map file: " + outFile + ".contigsMap.tsv");
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (Exception e) {
            throw new RuntimeException(e);
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

    private void validateScores() {
        // 0: scoreA - 100.00, 2: scoreB - (scoreA - 0.01), 1: scoreN - (scoreB - 0.01), -1: <= scoreN
        // if scoreB is 0 then scoreN must be 0
        if (scoreA < 0.0 || scoreA > 100.0) {
            Logger.error(CLASSNAME, "Score A must be between 0.0 and 100.0");
        }
        if (scoreB < 0.0 || scoreB > 100.0) {
            Logger.error(CLASSNAME, "Score B must be between 0.0 and 100.0");
        }
        if (scoreN < 0.0 || scoreN > 100.0) {
            Logger.error(CLASSNAME, "Score N must be between 0.0 and 100.0");
        }
        if (scoreA <= scoreB) {
            Logger.error(CLASSNAME, "Score A must be greater than Score B");
        }
        if (scoreB == scoreN){
            Logger.warning(CLASSNAME, "Score B is equal to Score N. There would be no alleles scored as het (1).");
            scoreN = scoreB;
        }
        if (scoreB == 0.0 && scoreN != 0.0) {
            Logger.warning(CLASSNAME, "Score B is not greater than Score N. There would be no alleles scored as missing (-1) or het (1).");
            scoreN = 0.0;
        }
    }
}
