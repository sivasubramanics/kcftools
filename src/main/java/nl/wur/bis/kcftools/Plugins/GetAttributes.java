package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.util.concurrent.Callable;

/***
 * This is a command line plugin that may extract attributes from the KCF files and write them to a TSV file
 * Attributes are:
 * 1. Number of observed kmers
 * 2. Number of variations
 * 3. Score
 * output file will be a TSV file with prefix.obs.tsv, prefix.var.tsv, prefix.score.tsv
 * The output files will have the following columns:
 * 1. window_id, sample1, sample2, ..., sampleN
 */
@Command(name = "getAttributes", description = "Extract attributes from KCF files")
public class GetAttributes implements Callable<Integer>, Runnable {
    @Option(names = {"-i", "--input"}, description = "KCF file name", required = true)
    private String kcfFile;
    @Option(names = {"-o", "--output"}, description = "Output file name prefix", required = true)
    private String outFile;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public GetAttributes() {
    }

    public GetAttributes(String kcfFile, String outFile) {
        this.kcfFile = kcfFile;
        this.outFile = outFile;
    }

    @Override
    public Integer call() throws Exception {
        getAttributes();
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /***
     * Main method to extract attributes from KCF file
     */
    private void getAttributes() {
        int[][] observedKmers;
        try (KCFReader reader = new KCFReader(kcfFile);
                BufferedWriter obsWriter = new BufferedWriter(new java.io.FileWriter(outFile + ".obs.tsv"));
                BufferedWriter varWriter = new BufferedWriter(new java.io.FileWriter(outFile + ".var.tsv"));
                BufferedWriter scoreWriter = new BufferedWriter(new java.io.FileWriter(outFile + ".score.tsv"));
             BufferedWriter totWriter = new BufferedWriter(new java.io.FileWriter(outFile + ".totalkmers.tsv"));
             BufferedWriter winlen = new BufferedWriter(new java.io.FileWriter(outFile + ".winlen.tsv"));
             BufferedWriter inDist = new BufferedWriter(new java.io.FileWriter(outFile + ".inDist.tsv"));
             BufferedWriter tailDist = new BufferedWriter(new java.io.FileWriter(outFile + ".tailDist.tsv"))) {
            Logger.info(CLASS_NAME, "Reading KCF file: " + kcfFile);
            KCFHeader header = reader.getHeader();
            String[] samples = header.getSamples();
            int nWindows = header.getWindowCount();
            obsWriter.write("window_id");
            varWriter.write("window_id");
            scoreWriter.write("window_id");
            inDist.write("window_id");
            tailDist.write("window_id");
            totWriter.write("window_id" + "\t" + "total_kmers");
            winlen.write("window_id" + "\t" + "window_length");
            for (String sample: samples){
                obsWriter.write("\t" + sample);
                varWriter.write("\t" + sample);
                scoreWriter.write("\t" + sample);
                inDist.write("\t" + sample);
                tailDist.write("\t" + sample);
            }
            obsWriter.newLine();
            varWriter.newLine();
            scoreWriter.newLine();
            totWriter.newLine();
            winlen.newLine();
            inDist.newLine();
            tailDist.newLine();
            for (Window window: reader){
                obsWriter.write(window.getWindowId());
                varWriter.write(window.getWindowId());
                scoreWriter.write(window.getWindowId());
                inDist.write(window.getWindowId());
                tailDist.write(window.getWindowId());
                for (String sample : samples) {
                    obsWriter.write("\t" + window.getObservedKmers(sample));
                    varWriter.write("\t" + window.getVariations(sample));
                    scoreWriter.write("\t" + String.format("%.2f", window.getScore(sample)));
                    inDist.write("\t" + window.getInnerDistance(sample));
                    tailDist.write("\t" + window.getTailDistance(sample));
                }
                totWriter.write(window.getWindowId() + "\t" + window.getTotalKmers());
                winlen.write(window.getWindowId() + "\t" + window.getEffLength());
                obsWriter.newLine();
                varWriter.newLine();
                scoreWriter.newLine();
                inDist.newLine();
                tailDist.newLine();
                totWriter.newLine();
                winlen.newLine();
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
//EOF