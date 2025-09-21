package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Callable;

/***
 * Convert the KCF windows to TSV file (just to validate with IBSpy results)
 * The output file will have the following columns:
 * 1. seqname
 * 2. start
 * 3. end
 * 4. total_kmers
 * 5. observed_kmers
 * 6. variations
 * 7. kmer_distance
 */
@Command(name = "kcf2tsv", description = "Convert KCF file to TSV file (IBSpy like)")
public class KCFToTSV implements Callable<Integer>, Runnable {
    private final String CLASSNAME = this.getClass().getSimpleName();
    @Option(names = {"-i", "--input"}, description = "KCF file name", required = true)
    private String kcfFile;
    @Option(names = {"-o", "--output"}, description = "Output file name Prefix", required = true)
    private String outFile;
    @Option(names = {"-s", "--sample"}, description = "Sample name", required = false)
    private String sampleName;

    public KCFToTSV() {
    }

    public KCFToTSV(String kcfFile, String outFile, String sampleName) {
        this.kcfFile = kcfFile;
        this.outFile = outFile;
        this.sampleName = sampleName;
    }

    @Override
    public Integer call() throws Exception {
        kcfToTSV();
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
     * Main method to convert KCF file to TSV file
     */
    private void kcfToTSV() {
        try (KCFReader reader = new KCFReader(kcfFile)) {
            KCFHeader header = reader.getHeader();
            Window[] windows = new Window[header.getWindowCount()];
            String[] querySamples;
            int winNum = 0;
            for (Window window : reader) {
                windows[winNum] = window;
                winNum++;
            }
            if (sampleName != null) {
                if (!header.hasSample(sampleName)) {
                    Logger.error(CLASSNAME, "Sample " + sampleName + " not found in KCF file");
                }
                querySamples = new String[]{sampleName};
            }
            else{
                querySamples = header.getSamples();
            }
            for (String sample : querySamples) {
                writeTSV(windows, sample);
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /***
     * Write the TSV file for the given sample
     */
    private void writeTSV(Window[] windows,String sample){
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile + "." + sample + ".tsv"))) {
            writer.write("window_id\tseqname\tstart\tend\teff_len\ttotal_kmers\tobserved_kmers\tvariations\tkmer_distance\tscore\n");
            for (Window window : windows) {
                if (window != null) {
                    writer.write(window.toTSV(sample) + "\n");
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
//EOF