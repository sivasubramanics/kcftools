package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.HelperFunctions;
import org.jetbrains.annotations.NotNull;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

/***
 * Class to read KCF files
 */
public class KCFReader implements Iterable<Window>, AutoCloseable {
    private final String filename;
    private KCFHeader header;
    private int windowId = 0;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public KCFReader(String filename) {
        this.filename = filename;
        this.header = null;
        HelperFunctions.log("info", CLASS_NAME, "Reading KCF file: " + filename);
    }

    public KCFHeader getHeader() {
        if (header == null) {
            try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
                StringBuilder headerBuilder = new StringBuilder();
                String line;
                while ((line = reader.readLine()) != null && line.startsWith("##")) {
                    headerBuilder.append(line).append("\n");
                }
                headerBuilder.append(line).append("\n");
                header = new KCFHeader(headerBuilder.toString());
            } catch (IOException e) {
                throw new RuntimeException("Error reading KCF file header", e);
            }
        }
        return header;
    }

    @Override
    public @NotNull Iterator<Window> iterator() {
        return new WindowIterator();
    }

    @Override
    public void close() throws Exception {
        // Nothing to do
    }

    private class WindowIterator implements Iterator<Window> {
        private final BufferedReader reader;
        private String nextLine;

        public WindowIterator() {
            try {
                reader = new BufferedReader(new FileReader(filename));
                // Skip header lines
                while ((nextLine = reader.readLine()) != null) {
                    if (nextLine.startsWith("#CHROM")) {
                        nextLine = reader.readLine();
                        break;
                    }
                }
            } catch (IOException e) {
                throw new RuntimeException("Error initializing KCF file reader", e);
            }
        }

        @Override
        public boolean hasNext() {
            return nextLine != null;
        }

        @Override
        public Window next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }

            try {
                String[] fields = nextLine.split("\t");
                String sequenceName = fields[0];
                int start = Integer.parseInt(fields[1]);
                int end = Integer.parseInt(fields[2]);
                int totalKmers = Integer.parseInt(fields[3]);

                Window window = new Window(windowId, sequenceName, start, end, totalKmers);
                windowId++;

                // Parse sample data
                for (int i = 0; i < getHeader().getSamples().length; i++) {
                    String[] sampleData = fields[i + 6].split(":");
                    window.addData(getHeader().getSamples()[i], Integer.parseInt(sampleData[2]), Integer.parseInt(sampleData[1]), Integer.parseInt(sampleData[3]), sampleData[0]);
                }

                nextLine = reader.readLine();
                return window;
            } catch (IOException e) {
                HelperFunctions.log("error", CLASS_NAME, "Error reading KCF file: " + filename);
                throw new RuntimeException("Error reading KCF file", e);
            }
        }
    }
}
// EOF