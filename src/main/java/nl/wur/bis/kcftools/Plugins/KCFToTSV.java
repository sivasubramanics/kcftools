package nl.wur.bis.kcftools.Plugins;


import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.Callable;

@Command(name = "kcf2tsv", description = "Convert KCF file to TSV file (IBSpy like)")
public class KCFToTSV implements Callable<Integer>, Runnable {
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
                    HelperFunctions.log("error", "KCFToTSV", "Sample " + sampleName + " not found in KCF file");
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

    private void writeTSV(Window[] windows,String sample){
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile + "." + sample + ".tsv"))) {
            writer.write("seqname\tstart\tend\ttotal_kmers\tobserved_kmers\tvariations\tkmer_distance\n");
            for (Window window : windows) {
                writer.write(window.toTSV(sample) + "\n");
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
