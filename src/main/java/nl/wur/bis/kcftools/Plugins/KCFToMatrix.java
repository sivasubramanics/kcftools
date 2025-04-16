package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
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

    @Option(names = {"-a", "--score_a"}, description = "lower score cut-off for reference allelse", required = false)
    private double scoreA = 95.0000;

    @Option(names = {"-b", "--score_b"}, description = "lower score cut-off for alternate allele", required = false)
    private double scoreB = 60.0000;

    @Option(names = {"-r", "--rdata"}, description = "Convert the matrix to RData", required = false)
    private boolean rdata = false;

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

    private void convertKCFToMatrix() {
        try (KCFReader reader = new KCFReader(inFile);
             BufferedWriter writer = new BufferedWriter(new java.io.FileWriter(outPrefix + ".matrix.tsv"));
             BufferedWriter mapWriter = new BufferedWriter(new java.io.FileWriter(outPrefix + ".map.tsv"))) {
            KCFHeader header = reader.getHeader();
            int[][] matrix = new int[header.getSamples().length][header.getWindowCount()];
            String[] samples = header.getSamples();
            int i = 0;
            for (Window window : reader) {
                mapWriter.write(window.getSequenceName() + "\t" + window.getStart() + "\t" + window.getEnd());
                mapWriter.newLine();
                for (int j = 0; j < samples.length; j++) {
                    // if score is between scoreCutOff and 100, then set the ibs to 0
                    if (window.getData().get(samples[j]).getScore() >= scoreA) {
                        matrix[j][i] = 0;
                    }
                    // if score is between 60 and scoreCutOff, then set the ibs to 2
                    else if (window.getData().get(samples[j]).getScore() >= scoreB && window.getData().get(samples[j]).getScore() < scoreA) {
                        matrix[j][i] = 2;
                    }
                    // if score is between 0 and 60, then set the ibs to 1
                    else if (window.getData().get(samples[j]).getScore() == 0) {
                        matrix[j][i] = 'N';
                    }
                    else {
                        matrix[j][i] = 1;
                    }
                }
                i++;
            }
            // write the matrix
            writer.write("sample");
            for (int k = 0; k < header.getWindowCount(); k++) {
                writer.write("\t" + k);
            }
            writer.newLine();
            for (int j = 0; j < samples.length; j++) {
                writer.write(samples[j]);
                for (int k = 0; k < header.getWindowCount(); k++) {
                    writer.write("\t" + matrix[j][k]);
                }
                writer.newLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (rdata) {
            convertGTmatrixToRdata(outPrefix + ".matrix.tsv", outPrefix + ".map.tsv");
        }
    }

    // convert matrix to Rdata
    private void convertGTmatrixToRdata(String matrixFile, String mapFile) {
        // check for Rscript is installed
        if (!HelperFunctions.isInstalled("Rscript")) {
            Logger.error(CLASSNAME, "Rscript is not installed. Please install Rscript and try again.");
        }
        // create random name for R script based on the current time point
        String rscriptName = "convertGTmatrixToRdata_" + System.currentTimeMillis() + ".R";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(rscriptName))) {
            writer.write("myGD <- read.table(\"" + matrixFile + "\", header = TRUE, sep = \"\\t\")\n");
            writer.write("save(myGD, file = \"" + matrixFile.replace(".tsv$", ".RData") + "\")\n");
            writer.write("myGM <- read.table(\"" + mapFile + "\", header = FALSE, sep = \"\\t\")\n");
            writer.write("save(myGM, file = \"" + mapFile.replace(".tsv$", ".RData") + "\")\n");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        // run the R script
        try {
            HelperFunctions.tryExec("Rscript " + rscriptName);
        }
        catch (Exception e) {
            Logger.error(CLASSNAME, "Error while running Rscript: " + e.getMessage());
        }
    }
}
// EOF