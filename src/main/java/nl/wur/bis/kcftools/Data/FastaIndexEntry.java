package nl.wur.bis.kcftools.Data;

import org.jetbrains.annotations.NotNull;

public class FastaIndexEntry implements Comparable<FastaIndexEntry> {
    private final String sequenceName;
    private final long offset;
    private final int length;
    private final int lineBases;
    private final int lineWidth;
    private final int seqId;

    public FastaIndexEntry(int seqId, String sequenceName, int length, long offset, int lineBases, int lineWidth) {
        this.seqId = seqId;
        this.sequenceName = sequenceName;
        this.offset = offset;
        this.length = length;
        this.lineBases = lineBases;
        this.lineWidth = lineWidth;
    }

    public String getSequenceName() {
        return sequenceName;
    }

    public long getOffset() {
        return offset;
    }

    public int getLength() {
        return length;
    }

    public int getLineBases() {
        return lineBases;
    }

    public int getLineWidth() {
        return lineWidth;
    }

    @Override
    public String toString() {
        return sequenceName + "\t" + length + "\t" + offset + "\t" + lineBases + "\t" + lineWidth;
    }

    @Override
    public int compareTo(@NotNull FastaIndexEntry o) {
        return Integer.compare(this.seqId, o.seqId);
    }

    public int getSeqId() {
        return seqId;
    }
}
//EOF