package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.KCFWriter;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.util.concurrent.Callable;

@Command(name = "splitKCF", description = "Split KCF file for each chromosome")
public class SplitKCF implements Callable<Integer>, Runnable {
    // in KCF file name
    @Option(names = {"-k", "--kcf"}, description = "KCF file name", required = true)
    private String kcfFile;
    // in output directory
    @Option(names = {"-o", "--output"}, description = "Output directory", required = true)
    private String outDir;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public SplitKCF() {
    }

    public SplitKCF(String kcfFile, String outDir) {
        this.kcfFile = kcfFile;
        this.outDir = outDir;
    }

    @Override
    public Integer call() throws Exception {
        splitKCF();
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

    private void splitKCF() {
        if (HelperFunctions.checkDirectoryExists(outDir)) {
            Logger.info(CLASS_NAME, "Output directory already exists: " + outDir);
        }
        else {
            Logger.info(CLASS_NAME, "Creating output directory: " + outDir);
            HelperFunctions.createDirectory(outDir);
        }
        try (KCFReader reader = new KCFReader(kcfFile)) {
            KCFHeader header = reader.getHeader();
            String[] samples = header.getSamples();
            String[] chromosomes = header.getContigs();
            for (String chromosome : chromosomes) {
                try (KCFWriter writer = new KCFWriter(outDir + "/" + chromosome + ".kcf")) {
                    writer.writeHeader(header);
                    for (Window window : reader) {
                        if (window.getSequenceName().equals(chromosome)) {
                            writer.writeWindow(window);
                        }
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error splitting KCF file", e);
        }
    }
}
