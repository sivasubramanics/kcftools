package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Logger;
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
        Logger.info(CLASS_NAME, "Reading KCF file:" + filename);
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

                Window window = new Window(fields, getHeader().getSamples(), getHeader().getWeights());

                nextLine = reader.readLine();
                return window;
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error reading KCF file: " + filename);
                throw new RuntimeException("Error reading KCF file", e);
            }
        }
    }
}
// EOF