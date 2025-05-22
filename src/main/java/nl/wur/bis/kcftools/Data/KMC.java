package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.BufferUnderflowException;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

import static nl.wur.bis.kcftools.Data.Kmer.binaryToKmer;

/***
 * This class represents a KMC object that is used to read a KMC database and query the database for kmer counts
 * This class is heavily based on the documentation here: https://github.com/refresh-bio/KMC/blob/master/API.pdf
 */
public class KMC implements AutoCloseable {
    private final String kmcPrefixFile;
    private final String kmcSuffixFile;
    private int kmerLength;
    private int mode;
    private int counterSize;
    private int lutPrefixLength;
    private int sufixLength;
    private int signatureLength;
    private int minCount;
    private int maxCount;
    private long totalKmers;
    private boolean bothStrands;
    private long singleLUTSize;
    private int version;
    // 120 mb
    int MAX_BYTE_COUNT = 1 << 24;
    private MappedByteBuffer[] suffixBuffers;
    private final Signature signatureReference;
    private final int record_size;
    private final long recordsPerPage;
    private int lutPrefixArraySize;

    private long[] prefixArray;
    private int[] signatureMap;
    private ThreadLocal<MappedByteBuffer[]> threadLocalBuffers;
    private byte[][] inMemorySuffixBuffers;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public KMC(String kmcDBName) throws IOException {
        this(kmcDBName, false);
    }

    public KMC(String kmcDBName, boolean inMemory) throws IOException {
        this.kmcPrefixFile = kmcDBName + ".kmc_pre";
        this.kmcSuffixFile = kmcDBName + ".kmc_suf";
        readPrefixFile(kmcPrefixFile);
        this.signatureReference = new Signature(signatureLength);
        record_size = counterSize + sufixLength / 4;
        recordsPerPage = MAX_BYTE_COUNT / record_size;
        if (inMemory) {
            preloadSuffixBuffers(kmcSuffixFile);
        }
        else {
            readSuffixBuffers(kmcSuffixFile);
            // add the thread-local buffers here
            threadLocalBuffers = ThreadLocal.withInitial(() -> {
                MappedByteBuffer[] localBuffers = new MappedByteBuffer[suffixBuffers.length];
                for (int i = 0; i < suffixBuffers.length; i++) {
                    localBuffers[i] = (MappedByteBuffer) suffixBuffers[i].duplicate();
                }
                return localBuffers;
            });
        }
        printSummary();
    }


    /***
     * Preload the suffix buffers into memory
     */
    private void preloadSuffixBuffers(String kmcSuffixFile) throws IOException {
        Logger.info(CLASS_NAME, "Loading KMC suffix file " + kmcSuffixFile + " into memory");
        Logger.warning(CLASS_NAME, "OutOfMemoryError may occur if the suffix file is too large");
        Logger.warning(CLASS_NAME, "Consider using the KCFTOOLS_HEAP_SIZE setting");
        int fullPageSize = MAX_BYTE_COUNT / record_size * record_size; // in bytes
        int numberOfPages = (int) ((totalKmers * record_size) / fullPageSize + ((totalKmers * record_size) % fullPageSize == 0 ? 0 : 1));

        inMemorySuffixBuffers = new byte[numberOfPages][];
        try (RandomAccessFile sufFile = new RandomAccessFile(kmcSuffixFile, "r")) {
            // first 4 bytes are the marker KMCS in the file
            sufFile.seek(4);
            for (int i = 0; i < numberOfPages; i++) {
                int pageSize = (int) (i == numberOfPages - 1 ? (totalKmers * record_size) % fullPageSize : fullPageSize);
                pageSize = pageSize == 0 ? fullPageSize : pageSize;
                inMemorySuffixBuffers[i] = new byte[pageSize];
                sufFile.readFully(inMemorySuffixBuffers[i]);
            }
        }
    }

    /***
     * Read the KMC prefix file
     */
    private void readPrefixFile(String kmcPrefixFile) throws IOException {
        Logger.info(CLASS_NAME, "Reading KMC prefix file " + kmcPrefixFile);
        try (RandomAccessFile raf = new RandomAccessFile(kmcPrefixFile, "r")) {
            long fileSize = raf.length();

            // memory map the file
            MappedByteBuffer buffer = raf.getChannel().map(FileChannel.MapMode.READ_ONLY, 0, fileSize);
            buffer.order(ByteOrder.LITTLE_ENDIAN);

            // read the header offset from the last 8 bytes
            buffer.position((int) (fileSize - 8));
            int headerOffset = buffer.getInt();

            // move to the header location in the file
            buffer.position((int) (fileSize - headerOffset - 8));

            // read header information
            kmerLength = buffer.getInt();
            mode = buffer.getInt();
            counterSize = buffer.getInt();
            lutPrefixLength = buffer.getInt();
            sufixLength = kmerLength - lutPrefixLength;
            signatureLength = buffer.getInt();
            minCount = buffer.getInt();
            maxCount = buffer.getInt();
            totalKmers = buffer.getLong();
            bothStrands = buffer.get() == 0;
            // skip uchar[3] padding
            buffer.position(buffer.position() + 3);
            // skip uint32[6] padding
            buffer.position(buffer.position() + 24);
            // get 4 bytes for version and check if its 0x200 or 0
            version = buffer.getInt();
            if (version != 0x200) {
                Logger.error(CLASS_NAME, "KMC version is not 0x200");
            }

            // go to position buffer.position((int) (fileSize - headerOffset - 8 - (signatureMapSize * 4)));
            long signatureMapSize = (1L << 2 * signatureLength) + 1;
            long signatureMapStart = fileSize - headerOffset - 8 - (signatureMapSize * 4);
            buffer.position((int) signatureMapStart);
            signatureMap = new int[(int) signatureMapSize];
            for (int i = 0; i < signatureMapSize; i++) {
                signatureMap[i] = buffer.getInt();
            }

            long prefixArrayStart = 4;
            lutPrefixArraySize = (1 << (2 * lutPrefixLength));
            singleLUTSize = lutPrefixArraySize * 8L;
            int numPrefixArrays = (int) ((signatureMapStart - 8 - 4) / (int) singleLUTSize);

            // go to prefixArrayStart
            buffer.position((int) prefixArrayStart);
            prefixArray = new long[numPrefixArrays * (1 << (2 * lutPrefixLength))];
            for (int i = 0; i < numPrefixArrays * (1 << (2 * lutPrefixLength)); i++) {
                prefixArray[i] = buffer.getLong();
            }

        } catch (IOException | BufferUnderflowException e) {
            Logger.error(CLASS_NAME, "Error reading prefix file " + kmcPrefixFile);
        }
    }

    /***
     * Read the suffix buffers from the kmc_suf file
     */
    private void readSuffixBuffers(String kmcSuffixFile){
        int full_page_size = MAX_BYTE_COUNT / record_size * record_size; // in bytes
        int number_of_pages = (int) ((totalKmers * record_size) / full_page_size + ((totalKmers * record_size) % full_page_size == 0 ? 0 : 1));
        Logger.info(CLASS_NAME, "MemoryMapping KMC suffix file " + kmcSuffixFile);
        try (RandomAccessFile suf_file = new RandomAccessFile(kmcSuffixFile, "r")) {
            // first 4 bytes are the marker KMCS in the file
            suf_file.seek(4);
            suffixBuffers = new MappedByteBuffer[number_of_pages];
            for (int i = 0; i < number_of_pages; ++i) {
                int page_size = (int) (i == number_of_pages - 1 ? (totalKmers * record_size) % full_page_size : full_page_size); // in bytes
                page_size = page_size == 0 ? full_page_size : page_size; // in bytes
                suffixBuffers[i] = suf_file.getChannel().map(FileChannel.MapMode.READ_ONLY, ((long)full_page_size) * i + 4, page_size);
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error reading suffix buffers from file " + kmcSuffixFile);
        }
    }

    /***
     * Dump the prefix array to a file [DEBUG function]
     */
    public void dumpPrefixArray(String prefixArrayFile) throws IOException {
        int columnSize = 1 << (2 * lutPrefixLength);
        Logger.info(CLASS_NAME, "Writing prefix array to file " + prefixArrayFile);
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(prefixArrayFile))) {
            for (int i = 0; i < prefixArray.length; i++) {
                writer.write(String.valueOf(prefixArray[i]));
                if (i % columnSize == columnSize - 1) {
                    writer.write("\n");
                } else {
                    writer.write("\t");
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing prefix array to file " + prefixArrayFile);
        }
    }

    /***
     * Dump the suffix buffers to a file [DEBUG function]
     */
    public void dumpSuffixBuffers(String suffixBufferFile) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(suffixBufferFile))) {
            // calculate the suffix size in bytes
            int suffixBytes = sufixLength / 4;
            // calculate entry size as suffix bytes plus counter size
            int entrySize = suffixBytes + counterSize;
            byte[] entry = new byte[entrySize];
            byte[] suffixByteArray;

            for (MappedByteBuffer mappedByteBuffer : suffixBuffers) {
                while (mappedByteBuffer.hasRemaining()) {
                    // get entry from buffer
                    mappedByteBuffer.get(entry);
                    int count = 0;
                    suffixByteArray = getSuffixFromEntry(entry);
                    count = getCountFromEntry(entry);

                    // convert suffix to string
                    String suffix = convertBytesToKmer(suffixByteArray, sufixLength);

                    // write to file
                    writer.write(suffix + " " + count);
                    writer.write("\n");
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing suffix buffer to file " + suffixBufferFile);
        }
    }

    /***
     * Dump the signature map to a file [DEBUG function]
     */
    public void dumpSignatureMap(String signatureMapFile) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(signatureMapFile))) {
            for (int i = 0; i < signatureMap.length; i++) {
                writer.write(i + " " + signatureMap[i] + "\n");
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing signature map to file " + signatureMapFile);
        }
    }


    /***
     * return the bothStrands flag
     */
    public boolean isBothStrands() {
        return bothStrands;
    }

    /***
     * Print a summary of the KMC object
     */
    public void printSummary() {
        Logger.info(CLASS_NAME, "==================== KMC INFO ====================");
        Logger.info(CLASS_NAME, String.format("%-25s: %s", "KMC prefix file", kmcPrefixFile));
        Logger.info(CLASS_NAME, String.format("%-25s: %s", "KMC suffix file", kmcSuffixFile));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Kmer length", kmerLength));
        Logger.info(CLASS_NAME, String.format("%-25s: %s", "Mode", mode));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Counter size", counterSize));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "LUT prefix length", lutPrefixLength));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Signature length", signatureLength));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Min count", minCount));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Max count", maxCount));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Total kmers", totalKmers));
        Logger.info(CLASS_NAME, String.format("%-25s: %b", "Both strands", bothStrands));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Signature map size", signatureMap.length));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Single LUT size", singleLUTSize));
        Logger.info(CLASS_NAME, String.format("%-25s: %d", "Prefix buffer size", prefixArray.length));
        Logger.info(CLASS_NAME, String.format("%-25s: %s", "Version", Integer.toHexString(version)));
        Logger.info(CLASS_NAME, "==================================================");

    }

    /***
     * Get the count of a kmer from the KMC database
     */
    public int getCount(Kmer kmer) {
        // get the signature, prefix, and suffix from the kmer
        long start, end;
        int signature = kmer.getSignature(this.signatureReference);
        int prefix = kmer.getPrefixFwd();
        byte[] suffix = kmer.getSuffixFwd();

        // calculate the start and end indices in the prefix array
        int signatureIndex = signatureMap[signature] * lutPrefixArraySize;
        start = prefixArray[signatureIndex + prefix];
        if (signatureIndex + prefix + 1 >= prefixArray.length) {
            end = totalKmers - 1;
        }
        else {
            end = prefixArray[signatureIndex + prefix + 1] - 1;
        }

        // perform binary search to find the matching suffix
        while (start <= end) {
            long mid = (start + end) / 2;
            byte[] entry = getEntry(mid); // retrieve the entry for the current index
            byte[] suffixFromEntry = getSuffixFromEntry(entry);

            int comparison = HelperFunctions.compareByteArray(suffix, suffixFromEntry);
            if (comparison < 0) {
                end = mid - 1; // move the search range down
            } else if (comparison > 0) {
                start = mid + 1; // move the search range up
            } else {
                return getCountFromEntry(entry); // Match found, return count
            }
        }

        return 0;
    }

    /***
     * Checks if the given kmer exists in the KMC database.
     */
    public boolean isExist(Kmer kmer){
        return getCount(kmer) != 0;
    }


    /***
     * Converts a byte array to a string.
     */
    private String convertBytesToKmer(byte[] suffixBytes, int suffixLength) {
        StringBuilder suffix = new StringBuilder(suffixLength);
        for (int i = 0; i < suffixLength; ++i) {
            int baseIndex = (suffixBytes[i / 4] >> (6 - 2 * (i % 4))) & 0x03;
            suffix.append(baseToChar(baseIndex));
        }
        return suffix.toString();
    }

    /***
     * Converts a nucleotide base index to its character representation.
     * @param baseIndex Base index (0, 1, 2, 3 for A, C, G, T)
     * @return Nucleotide character
     */
    private char baseToChar(int baseIndex) {
        return switch (baseIndex) {
            case 0b00 -> 'A';
            case 0b01 -> 'C';
            case 0b10 -> 'G';
            case 0b11 -> 'T';
            default -> throw new IllegalArgumentException("Invalid base index: " + baseIndex);
        };
    }

    /***
     * Get the suffix entry from the kmc_suf file (THREAD SAFE, so it will be a bit slower (very minimal))
     */
    public byte[] getEntry(long index) {
        long page = index / recordsPerPage;
        int offset = (int) (index % recordsPerPage);
        int recordOffset = offset * record_size;
        byte[] entry = new byte[record_size];

        if (inMemorySuffixBuffers != null) {
            System.arraycopy(inMemorySuffixBuffers[(int) page], recordOffset, entry, 0, record_size);
        }
        else{
            MappedByteBuffer localBuffer = threadLocalBuffers.get()[(int) page];
            localBuffer.position(recordOffset);
            localBuffer.get(entry);
        }
        return entry;
    }

    /***
     * Get the suffix from the suffix entry from the kmc_suf file
     */
    public byte[] getSuffixFromEntry(byte[] entry) {
        byte[] suffix = new byte[sufixLength / 4];
        System.arraycopy(entry, 0, suffix, 0, sufixLength / 4);
        return suffix;
    }

    /***
     * Get the count from the suffix entry from the kmc_suf file
     */
    public int getCountFromEntry(byte[] entry) {
        int count = 0;
        for (int i = 0; i < counterSize; i++) {
            count |= (entry[sufixLength / 4 + i] & 0xFF) << (i * 8);
        }
        return count;
    }

    public int getKmerLength() {
        return kmerLength;
    }


    public int getPrefixLength() {
        return lutPrefixLength;
    }

    @Override
    public void close() {
        if (suffixBuffers != null) {
            Arrays.fill(suffixBuffers, null);
        }
        suffixBuffers = null;
        prefixArray = null;
        signatureMap = null;
    }

    /***
     * dump the kmer table to a file (handle with care as this could create huge file) [DEBUG function]
     * @param kmerTableFile
     * @throws IOException
     */
    public void dumpKmerTable(String kmerTableFile) throws IOException {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(kmerTableFile))) {
            // iterate the whole prefixArray and print the prefix and suffix and count for each kmer
            for (int i = 0; i < prefixArray.length; i++) {
                long start = prefixArray[i];
                long end = (i == prefixArray.length - 1) ? totalKmers - 1 : prefixArray[i + 1] - 1;
                // get the prefix from the i value where i % singleLUTSize is the prefix int
                int prefix = i % lutPrefixArraySize;
                // conver int to string
                char[] pre = binaryToKmer(prefix, lutPrefixLength);
                String prefixString = new String(pre);
                for (long j = start; j <= end; j++) {
                    byte[] entry = getEntry(j);
                    byte[] suffix = getSuffixFromEntry(entry);
                    int count = getCountFromEntry(entry);
                    String suf = convertBytesToKmer(suffix, sufixLength);
                    writer.write( prefixString + suf + "\t" + count + "\n");
                }
            }

        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing kmer table to file " + kmerTableFile);
        }
    }
}
//EOF