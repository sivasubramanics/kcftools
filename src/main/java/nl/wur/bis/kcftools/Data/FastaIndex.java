package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;

import java.io.*;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/***
 * This class represents a FastaIndex object which is used to index a fasta file (somewhere similar to samtools faidx)
 */
public class FastaIndex implements AutoCloseable {
    private final Map<String, FastaIndexEntry> index;
    private final List<String> sequenceNames;
    private final MappedByteBuffer[] buffer;
    private final Object bufferLock = new Object();
    private final String CLASS_NAME = this.getClass().getSimpleName();

    public FastaIndex(String fastaFilePath) throws IOException {
        File fastaFile = new File(fastaFilePath);
        String faiFilePath = fastaFilePath + ".faidx";
        File indexFile = new File(faiFilePath);

        if (!indexFile.exists() || HelperFunctions.isOlder(indexFile, fastaFile)) {
            Logger.info(CLASS_NAME, "Generating/Updating index file: " + faiFilePath);
            generateIndexFile(fastaFile, indexFile);
        } else {
            Logger.info(CLASS_NAME, "Using existing index file: " + faiFilePath);
        }

        // sort this.index by entry.seqId
        this.index = sortByValue(Collections.unmodifiableMap(loadIndex(indexFile)));
        // Create an immutable list of sequence names in the same order as in the index
        this.sequenceNames = Collections.unmodifiableList(new ArrayList<>(index.keySet()));

        try (RandomAccessFile randomAccessFile = new RandomAccessFile(fastaFile, "r");
             FileChannel fileChannel = randomAccessFile.getChannel()) {

            long fileSize = fileChannel.size();
            int numBuffers = this.sequenceNames.size();
            this.buffer = new MappedByteBuffer[numBuffers];

            for (int i = 0; i < numBuffers; i++) {
                String sequenceName = this.sequenceNames.get(i);
                FastaIndexEntry entry = this.index.get(sequenceName);

                long sequenceStartOffset = entry.getOffset();
                long sequenceEndOffset;
                if (i + 1 < numBuffers) {
                    FastaIndexEntry nextEntry = this.index.get(this.sequenceNames.get(i + 1));
                    sequenceEndOffset = nextEntry.getOffset();
                } else {
                    sequenceEndOffset = fileSize;
                }

                // Offset length of the sequence
                int sequenceOffsetLength = (int) (sequenceEndOffset - sequenceStartOffset);

                // Memory-map for each sequence
                if (sequenceOffsetLength > 0) {
                    this.buffer[i] = fileChannel.map(FileChannel.MapMode.READ_ONLY, sequenceStartOffset, sequenceOffsetLength);
                } else {
                    this.buffer[i] = null;
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error memory-mapping fasta file: " + e.getMessage());
            throw e;
        }
    }

    /***
     * Load the index file into a map
     */
    private Map<String, FastaIndexEntry> loadIndex(File indexFile) throws IOException {
        Map<String, FastaIndexEntry> tempIndex = new ConcurrentHashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(indexFile))) {
            String line;
            int seqId = 0;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                FastaIndexEntry entry = new FastaIndexEntry(seqId,
                        fields[0],
                        Integer.parseInt(fields[1]),
                        Long.parseLong(fields[2]),
                        Integer.parseInt(fields[3]),
                        Integer.parseInt(fields[4])
                );
                seqId++;
                if (tempIndex.putIfAbsent(entry.getSequenceName(), entry) != null) {
                    Logger.error(CLASS_NAME, "Duplicate sequence name in index: " + entry.getSequenceName());
                }
            }
        }
        return tempIndex;
    }

    /***
     * Get a FastaIndexEntry by sequence name
     */
    public FastaIndexEntry getEntry(String sequenceName) {
        return index.get(sequenceName);
    }

    /***
     * Get a sequence from the fasta file by name
     */
    public String getSequence(String sequenceName) {
        return getSequence(sequenceName, 0, getEntry(sequenceName).getLength());
    }

    /***
     * Get a subsequence from the fasta file by name, start and length
     */
    public String getSequence(String name, int start, int length) {
        int end = start + length;
        FastaIndexEntry entry = getEntry(name);
        if (entry == null) {
            Logger.error(CLASS_NAME, "Sequence not found in index: " + name);
            return null;
        }

        // Validate start and end
        if (start < 0 || end > entry.getLength() || start >= end) {
            Logger.error(CLASS_NAME, "Invalid range: " + start + "-" + end + " for sequence: " + name);
            return null;
        }

        StringBuilder sequence = new StringBuilder(end - start);

        synchronized (bufferLock) {
            try {

                // get the buffer for the sequence, and calculate the start and end positions in the buffer
                MappedByteBuffer buf = buffer[entry.getSeqId()];

                int lineBases = entry.getLineBases();
                int lineWidth = entry.getLineWidth();

                // get the line at which the start position is located
                int startLine = start / lineBases;
                int startLineBaseIndex = start % lineBases;

                // get the offset of the start position in the line
                long startOffset = ((long) startLine * lineWidth) + startLineBaseIndex;

                buf.position((int) (startOffset));

                // read the sequence from the buffer
                int basesToExtract = end - start;
                while (basesToExtract > 0) {
                    // get number of bases to read in the current line (or remaining bases in the line)
                    int remainingBasesInLine = lineBases - (startLineBaseIndex % lineBases);
                    int basesToRead = Math.min(basesToExtract, remainingBasesInLine);

                    // read the bases
                    for (int i = 0; i < basesToRead; i++) {
                        sequence.append((char) buf.get());
                    }

                    // move to the next line
                    buf.position(buf.position() + (lineWidth - lineBases));
                    basesToExtract -= basesToRead;

                    // reset the startLineBaseIndex
                    startLineBaseIndex = 0;
                }
            } catch (IllegalArgumentException | IndexOutOfBoundsException e) {
                Logger.error(CLASS_NAME, "Error reading sequence: " + e.getMessage());
                return null;
            }
        }

        return sequence.toString();
    }

    /***
     * Check if the sequence name is present in the index
     */
    public boolean containsSequence(String sequenceName) {
        return index.containsKey(sequenceName);
    }

    @Override
    public void close() throws IOException {
        for (MappedByteBuffer buf : buffer) {
            if (buf != null) {
                buf.clear();
            }
        }
    }

    /***
     * Get a Fasta object for a given sequence name
     */
    public Fasta getFasta(String tr) {
        synchronized (bufferLock) {
            int seqId = sequenceNames.indexOf(tr);
            return new Fasta(seqId, tr, getSequence(tr));
        }
    }

    /***
     * Get the number of sequences in the fasta file
     */
    public int length() {
        return sequenceNames.size();
    }

    /***
     * Get the sequence names
     */
    public List<String> getSequenceNames() {
        return sequenceNames;
    }

    /***
     * Get the length of a sequence
     */
    public int getSequenceLength(String name) {
        FastaIndexEntry entry = getEntry(name);
        if (entry == null) {
            Logger.error(CLASS_NAME, "Sequence not found in index: " + name);
            return -1;
        }
        return entry.getLength();
    }

    /***
     * Generate an index file for the fasta file
     */
    private void generateIndexFile(File fastaFile, File indexFile) throws IOException {
        if (HelperFunctions.isCompressed(fastaFile)) {
            Logger.error(CLASS_NAME, "Fasta file is compressed. Please decompress before indexing: " + fastaFile);
            throw new IllegalArgumentException("Fasta file is compressed. Please decompress before indexing: " + fastaFile);
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile));
             BufferedWriter writer = new BufferedWriter(new FileWriter(indexFile))) {

            String line;
            long offset = 0;
            String currentName = null;
            int lineBases = 0;
            int lineWidth = 0;
            int seqLength = 0;
            long sequenceStartOffset = 0;
            boolean inSequence = false;

            while ((line = reader.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (inSequence) {
                        writer.write(currentName + "\t" + seqLength + "\t" + sequenceStartOffset + "\t" + lineBases + "\t" + lineWidth + "\n");
                    }
                    currentName = line.substring(1).split(" ")[0];
                    offset += line.length() + System.lineSeparator().length();
                    sequenceStartOffset = offset;
                    inSequence = true;
                    seqLength = 0;
                } else {
                    if (seqLength == 0) {
                        lineBases = line.length();
                        lineWidth = line.length() + System.lineSeparator().length();
                    }
                    seqLength += line.length();
                    offset += line.length() + System.lineSeparator().length();
                }
            }
            if (inSequence && currentName != null) {
                writer.write(currentName + "\t" + seqLength + "\t" + sequenceStartOffset + "\t" + lineBases + "\t" + lineWidth + "\n");
            }
        }
    }

    /***
     * Sort the index Map by value
     */
    public static Map<String, FastaIndexEntry> sortByValue(Map<String, FastaIndexEntry> map) {
        List<Map.Entry<String, FastaIndexEntry>> entryList = new ArrayList<>(map.entrySet());
        entryList.sort(Map.Entry.comparingByValue());
        Map<String, FastaIndexEntry> sortedMap = new LinkedHashMap<>();
        for (Map.Entry<String, FastaIndexEntry> entry : entryList) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }
        return sortedMap;
    }
}
//EOF