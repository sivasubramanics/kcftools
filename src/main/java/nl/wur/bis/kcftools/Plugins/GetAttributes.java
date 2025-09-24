package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.Callable;

/***
 * Command-line plugin to extract attributes from KCF files and write them to TSV files.
 * Supported attributes:
 *   - obs        : observed kmers
 *   - var        : variations
 *   - kd         : mean kmer count
 *   - score      : score
 *   - totalkmers : total kmers per window
 *   - winlen     : effective window length
 *   - inDist     : inner distance
 *   - tailDist   : tail distance
 *sdfasdd
 * Example:
 *   kcftools getAttributes -i input.kcf -o outPrefix -a obs var score
 *   kcftools getAttributes -i input.kcf -o outPrefix   (extracts all attributes)
 */
@Command(name = "getAttributes", description = "Extract attributes from KCF files")
public class GetAttributes implements Callable<Integer>, Runnable {

    @Option(names = {"-i", "--input"}, description = "KCF file name", required = true)
    private String kcfFile;

    @Option(names = {"-o", "--output"}, description = "Output file name prefix", required = true)
    private String outFile;

    @Option(
            names = {"-a", "--attributes"},
            split = ",",
            description = {
                    "Attributes to extract. Default: all",
                    "  - obs        : observed kmers",
                    "  - var        : variations",
                    "  - kd         : mean kmer count",
                    "  - score      : score",
                    "  - totalkmers : total kmers per window",
                    "  - winlen     : effective window length",
                    "  - inDist     : inner distance",
                    "  - tailDist   : tail distance"
            }
    )
    private List<String> attributes;


    private final String CLASS_NAME = this.getClass().getSimpleName();

    // Define all available attributes once
    private static final List<String> ALL_ATTRIBUTES = Arrays.asList(
            "obs", "var", "kd", "score", "totalkmers", "winlen", "inDist", "tailDist"
    );

    @Override
    public Integer call() throws Exception {
        extractAttributes();
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
     * Extract selected attributes and write to TSV.
     */
    private void extractAttributes() {
        try (KCFReader reader = new KCFReader(kcfFile)) {

            KCFHeader header = reader.getHeader();
            String[] samples = header.getSamples();

            // Default to all attributes if none were provided
            List<String> attrsToExtract = (attributes == null || attributes.isEmpty())
                    ? ALL_ATTRIBUTES
                    : attributes;

            Logger.info(CLASS_NAME, "Extracting attributes: " + String.join(", ", attrsToExtract));

            // Prepare writers only for requested attributes
            Map<String, BufferedWriter> writers = new HashMap<>();
            for (String attr : attrsToExtract) {
                if (!ALL_ATTRIBUTES.contains(attr)) {
                    throw new IllegalArgumentException("Unsupported attribute: " + attr);
                }
                String filename = outFile + "." + attr + ".tsv";
                writers.put(attr, new BufferedWriter(new FileWriter(filename)));
            }

            // Write headers
            for (Map.Entry<String, BufferedWriter> entry : writers.entrySet()) {
                String attr = entry.getKey();
                BufferedWriter w = entry.getValue();

                if (attr.equals("totalkmers")) {
                    w.write("window_id\ttotal_kmers");
                } else if (attr.equals("winlen")) {
                    w.write("window_id\twindow_length");
                } else {
                    w.write("window_id");
                    for (String sample : samples) {
                        w.write("\t" + sample);
                    }
                }
                w.newLine();
            }

            // Write data per window
            for (Window window : reader) {
                for (Map.Entry<String, BufferedWriter> entry : writers.entrySet()) {
                    String attr = entry.getKey();
                    BufferedWriter w = entry.getValue();

                    switch (attr) {
                        case "obs":
                            writeSampleValues(w, window.getWindowId(), samples, window::getObservedKmers);
                            break;
                        case "var":
                            writeSampleValues(w, window.getWindowId(), samples, window::getVariations);
                            break;
                        case "kd":
                            writeSampleValues(w, window.getWindowId(), samples,
                                    s -> String.format("%.2f", window.getMeanKmerCount(s)));
                            break;
                        case "score":
                            writeSampleValues(w, window.getWindowId(), samples,
                                    s -> String.format("%.2f", window.getScore(s)));
                            break;
                        case "inDist":
                            writeSampleValues(w, window.getWindowId(), samples, window::getInnerDistance);
                            break;
                        case "tailDist":
                            writeSampleValues(w, window.getWindowId(), samples, window::getTailDistance);
                            break;
                        case "totalkmers":
                            w.write(window.getWindowId() + "\t" + window.getTotalKmers());
                            w.newLine();
                            break;
                        case "winlen":
                            w.write(window.getWindowId() + "\t" + window.getEffLength());
                            w.newLine();
                            break;
                    }
                }
            }

            // Close writers
            for (BufferedWriter w : writers.values()) {
                w.close();
            }

        } catch (Exception e) {
            throw new RuntimeException("Error extracting attributes", e);
        }
    }

    /***
     * Helper to write values for each sample.
     */
    private void writeSampleValues(BufferedWriter w, String windowId,
                                   String[] samples, java.util.function.Function<String, Object> getter) throws Exception {
        w.write(windowId);
        for (String sample : samples) {
            w.write("\t" + getter.apply(sample));
        }
        w.newLine();
    }
}
// EOF