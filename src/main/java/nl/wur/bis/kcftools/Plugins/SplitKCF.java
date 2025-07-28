package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.KCFWriter;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.*;

/***
 * This is a command line plugin that splits the KCF file into separate files for each chromosome
 */
@Command(name = "splitKCF", description = "Split KCF file for each chromosome")
public class SplitKCF implements Callable<Integer>, Runnable {
    // in KCF file name
    @Option(names = {"-k", "--kcf"}, description = "KCF file name", required = true)
    private String kcfFile;
    // in output directory
    @Option(names = {"-o", "--output"}, description = "Output directory", required = true)
    private String outDir;
    // number of threads to use
    @Option(names = {"-t", "--threads"}, description = "Number of threads to use (default: auto-detected)", defaultValue = "2")
    private int nThreads = 2;

    private static final int MAX_OPEN_WRITERS = 100;

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
        } else {
            Logger.info(CLASS_NAME, "Creating output directory: " + outDir);
            HelperFunctions.createDirectory(outDir);
        }

        LinkedHashMap<String, KCFWriter> writerCache = new LinkedHashMap<>(MAX_OPEN_WRITERS, 0.75f, true) {
            @Override
            protected boolean removeEldestEntry(Map.Entry<String, KCFWriter> eldest) {
                if (size() > MAX_OPEN_WRITERS) {
                    eldest.getValue().close();
                    return true;
                }
                return false;
            }
        };

        try (KCFReader reader = new KCFReader(kcfFile)) {
            KCFHeader header = reader.getHeader();

            for (Window window : reader) {
                String chromosome = window.getSequenceName();
                KCFWriter writer = writerCache.get(chromosome);

                if (writer == null) {
                    writer = new KCFWriter(outDir + "/" + chromosome + ".kcf");
                    writer.writeHeader(header);
                    writerCache.put(chromosome, writer);
                }

                writer.writeWindow(window);
            }
        } catch (Exception e) {
            throw new RuntimeException("Error splitting KCF file", e);
        } finally {
            for (KCFWriter writer : writerCache.values()) {
                writer.close();
            }
        }
    }
}
//EOF