package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
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
                    // matrix[sample][window]
                    matrix[j][i] = window.getIbs(samples[j]);
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
    }
}
// EOF