package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Logger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/***
 * Class to write KCF files
 */
public class KCFWriter implements AutoCloseable {
    private KCFHeader header;
    private final BufferedWriter writer;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public KCFWriter(String filename) throws IOException {
        this.header = null;
        this.writer = new BufferedWriter(new FileWriter(filename));
        Logger.info(CLASS_NAME, "Writing KCF file: " + filename);
    }

    public void writeKCF(KCFHeader header, Iterable<Window> windows) throws IOException {
        if (header.getSamples() == null) {
            throw new IllegalArgumentException("Sample name not set in KCF header");
        }
        writeHeader(header);
        writeWindows(windows);
    }

    public void writeHeader(KCFHeader header) {
        this.header = header;
        try {
            writer.write(header.toString());
        } catch (IOException e) {
            throw new RuntimeException("Error writing KCF file header", e);
        }
    }

    public void writeWindow(Window window) {
        try {
            if (header.getSamples().length != window.getSamples().length) {
                throw new IllegalArgumentException("Number of samples in header does not match number of values in window");
            }
            writer.write(window.toString());
            writer.newLine();
        } catch (IOException e) {
            throw new RuntimeException("Error writing KCF file window", e);
        }
    }

    public void writeWindows(Iterable<Window> windows) throws IOException {
        for (Window window : windows) {
            writeWindow(window);
        }
        writer.flush();
    }

    public void close() {
        try {
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException("Error closing KCF file", e);
        }
    }
}
