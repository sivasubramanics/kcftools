package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.Arrays;

public class Kmer implements Comparable<Kmer> {
    private final int kmerLength;
    private final int prefixLength;
    private final int suffixLength;
    private final boolean isCanonical;
    private long[] kmerLong;

    private Integer prefixFwd;
    private byte[] suffixFwd;
    private Integer prefixRev;
    private byte[] suffixRev;
    private String reverseComplementKmer;

    private static final String CLASS_NAME = Kmer.class.getName();

    public Kmer(String kmer) {
        this(kmer.toCharArray(), 5, true);
    }

    public Kmer(String kmer, int prefixLength) {
        this(kmer.toCharArray(), prefixLength, true);
    }

    public Kmer(String kmer, int prefixLength, boolean isCanonical) {
        this(kmer.toCharArray(), prefixLength, isCanonical);
    }

    public Kmer(char @NotNull [] kmer, int prefixLength, boolean isCanonical) {
        this.kmerLength = kmer.length;
        this.prefixLength = prefixLength;
        this.suffixLength = this.kmerLength - prefixLength;
        this.isCanonical = isCanonical;
        this.kmerLong = kmerToLong(kmer);
        if (isCanonical) {
            // get the lowest of the kmer and its reverse complement
            long[] revComp = getReverseComplement(this.kmerLong, this.kmerLength);
            if (compareLongArrays(this.kmerLong, revComp) > 0) {
                this.kmerLong = revComp;
            }
        }
    }

    public Kmer(long[] kmerLong, int kmerLength, int prefixLength, boolean isCanonical) {
        this.kmerLength = kmerLength;
        this.prefixLength = prefixLength;
        this.suffixLength = this.kmerLength - prefixLength;
        this.isCanonical = isCanonical;
        this.kmerLong = Arrays.copyOf(kmerLong, kmerLong.length);
        if (isCanonical) {
            // get the lowest of the kmer and its reverse complement
            long[] revComp = getReverseComplement(this.kmerLong, this.kmerLength);
            if (compareLongArrays(this.kmerLong, revComp) > 0) {
                this.kmerLong = revComp;
            }
        }
    }

    public Kmer(Kmer kmer, boolean isCanonical) {
        // copy the kmer and set the canonical status
        this.kmerLength = kmer.getKmerLength();
        this.prefixLength = kmer.getPrefixLength();
        this.suffixLength = kmer.getSuffixLength();
        this.isCanonical = isCanonical;
        this.kmerLong = Arrays.copyOf(kmer.kmerLong, kmer.kmerLong.length);
        if (isCanonical) {
            // get the lowest of the kmer and its reverse complement
            long[] revComp = getReverseComplement(this.kmerLong, this.kmerLength);
            if (compareLongArrays(this.kmerLong, revComp) > 0) {
                this.kmerLong = revComp;
            }
        }
    }

    public int getKmerLength() {
        return this.kmerLength;
    }

    public int getPrefixLength() {
        return this.prefixLength;
    }

    public int getSuffixLength() {
        return this.suffixLength;
    }

    public String getSignatureString(Signature signature) {
        return longToKmer(this.getSignature(signature), signature.getSignLength());
    }

    public long[] getKmerLong() {
        return this.kmerLong;
    }

    /***
     * This method is used to get the signature of a kmer
     */
    public int getSignature(Signature signature) {
        int signatureLength = signature.getSignLength();
        int currentSignature = extractIntFromBits(this.kmerLong, 0, signatureLength);
        int minSignature = signature.getSignature(currentSignature);

        for (int i = 1; i <= this.kmerLength - signatureLength; i++) {
            currentSignature = ((currentSignature << 2) & ((1 << (2 * signatureLength)) - 1)) |
                    extractIntFromBits(this.kmerLong, i + signatureLength - 1, 1);
            if (signature.getSignature(currentSignature) < minSignature) {
                minSignature = signature.getSignature(currentSignature);
            }
        }
        return minSignature;
    }

    /***
     * Lazily calculates and returns the forward prefix as an integer.
     */
    public int getPrefixFwd() {
        if (prefixFwd == null) {
            prefixFwd = extractIntFromBits(kmerLong, 0, prefixLength);
        }
        return prefixFwd;
    }

    /***
     * Lazily calculates and returns the forward suffix as a byte array.
     */
    public byte[] getSuffixFwd() {
        if (suffixFwd == null) {
            suffixFwd = extractSuffix(kmerLong, kmerLength, prefixLength);
        }
        return suffixFwd;
    }

    /***
     * extract suffix from the kmer long and convert that to a byte array.
     */
    private byte[] extractSuffix(long[] kmerLong, int kmerLength, int prefixLength) {
        int suffixLength = kmerLength - prefixLength;
        byte[] suffix = new byte[(suffixLength + 3) / 4];

        int bitStart = 2 * prefixLength;
        int totalBits = 2 * kmerLength;

        for (int i = 0; i < suffixLength; ++i) {
            int bitIndex = bitStart + i * 2;

            int wordIndex = bitIndex / 64;
            int bitOffset = bitIndex % 64;

            int baseBits;
            if (bitOffset <= 62) {
                baseBits = (int) ((kmerLong[wordIndex] >> (62 - bitOffset)) & 0x3);
            } else {
                // Bits are split between two longs
                int highBits = (int) ((kmerLong[wordIndex] & 1L) << 1);
                int lowBits = (int) ((kmerLong[wordIndex + 1] >> 63) & 1L);
                baseBits = highBits | lowBits;
            }

            suffix[i / 4] |= (byte) (baseBits << ((3 - i % 4) * 2));
        }

        return suffix;
    }

    /***
     * extract the prefix from the kmer long and return that as an integer. [miscellaneous]
     * @return
     */
    public int getPrefixRev() {
        if (prefixRev == null) {
            long[] rev = getReverseComplement(kmerLong, kmerLength);
            prefixRev = extractIntFromBits(rev, 0, prefixLength);
        }
        return prefixRev;
    }

    /***
     * extract the suffix from the kmer long and return that as a byte array. [miscellaneous]
     * @return
     */
    public byte[] getSuffixRev() {
        if (suffixRev == null) {
            long[] rev = getReverseComplement(kmerLong, kmerLength);
            suffixRev = Arrays.copyOfRange(kmerToBytes(rev, kmerLength), prefixLength, kmerLength);
        }
        return suffixRev;
    }

    /***
     * get the reverse complement of the kmer long and return that as a string. [miscellaneous]
     * @return
     */
    public String getReverseComplementKmer() {
        if (reverseComplementKmer == null) {
            reverseComplementKmer = longToKmer(getReverseComplement(kmerLong, kmerLength), kmerLength);
        }
        return reverseComplementKmer;
    }

    /***
     * extract the integer value from the kmer long. used in the getSignature method to extract the best signature from the long kmer
     * @param bitArray
     * @param startBase
     * @param lengthBases
     * @return
     */
    private int extractIntFromBits(long[] bitArray, int startBase, int lengthBases) {
        int result = 0;
        for (int i = 0; i < lengthBases; i++) {
            int bitIndex = (startBase + i) * 2;
            int arrayIndex = bitIndex / 64;
            int offset = bitIndex % 64;
            int baseBits;
            if (offset <= 62) {
                baseBits = (int) ((bitArray[arrayIndex] >>> (62 - offset)) & 0x3);
            } else {
                int remaining = 64 - offset;
                long firstPart = (bitArray[arrayIndex] & ((1L << remaining) - 1)) << (2 - remaining);
                long secondPart = (bitArray[arrayIndex + 1] >>> (64 - (2 - remaining)));
                baseBits = (int) (firstPart | secondPart) & 0x3;
            }
            result = (result << 2) | baseBits;
        }
        return result;
    }

    /***
     * convert the kmer to a long array. used in the constructor to convert the kmer char array to a long array.
     * @param kmer
     * @return
     */
    public static long[] kmerToLong(char[] kmer) {
        int totalBits = kmer.length * 2;
        int arraySize = (totalBits + 63) / 64;
        long[] result = new long[arraySize];

        for (int i = 0; i < kmer.length; i++) {
            int val = baseToBits(kmer[i]);
            int bitIndex = i * 2;
            int arrayIndex = bitIndex / 64;
            int bitOffset = bitIndex % 64;

            result[arrayIndex] |= ((long) val) << (62 - bitOffset);

            // If split between two longs
            if (bitOffset > 62) {
                result[arrayIndex + 1] |= ((long) val) >>> (bitOffset - 62);
            }
        }

        return result;
    }

    /***
     * convert the kmer long array to a byte array. used in the getSuffixFwd and getSuffixRev methods to convert the kmer long array to a byte array.
     * @param kmerLong
     * @param kmerLength
     * @return
     */
    public static byte[] kmerToBytes(long[] kmerLong, int kmerLength) {
        byte[] bases = new byte[kmerLength];

        for (int i = 0; i < kmerLength; i++) {
            int bitIndex = i * 2;
            int arrayIndex = bitIndex / 64;
            int bitOffset = bitIndex % 64;

            int bits;
            if (bitOffset <= 62) {
                bits = (int) ((kmerLong[arrayIndex] >>> (62 - bitOffset)) & 0b11);
            } else {
                int bitsInFirstLong = 64 - bitOffset;
                long part1 = (kmerLong[arrayIndex] & ((1L << bitsInFirstLong) - 1)) << (2 - bitsInFirstLong);
                long part2 = (kmerLong[arrayIndex + 1] >>> (64 - (2 - bitsInFirstLong)));
                bits = (int) ((part1 | part2) & 0b11);
            }

            bases[i] = (byte) bits;
        }

        return bases;
    }

    /***
     * convert the base to bits. used in the kmerToLong method to convert the kmer char array to a long array.
     * @param base
     * @return
     */
    private static int baseToBits(char base) {
        return switch (base) {
            case 'A' -> 0b00;
            case 'C' -> 0b01;
            case 'G' -> 0b10;
            case 'T' -> 0b11;
            default -> throw new IllegalArgumentException("Invalid nucleotide in k-mer: " + base);
        };
    }

    /***
     * get the reverse complement of the kmer long array. used in the constructor to get the reverse complement of the kmer long array.
     * @param binaryKmer
     * @param kmerLength
     * @return
     */
    public static long[] getReverseComplement(long[] binaryKmer, int kmerLength) {
        int totalBits = kmerLength * 2;
        long[] reverse = new long[binaryKmer.length];

        for (int i = 0; i < kmerLength; i++) {
            // get the i-th base (forward)
            int bitIndex = i * 2;
            int arrayIndex = bitIndex / 64;
            int bitOffset = bitIndex % 64;

            int bits;
            if (bitOffset <= 62) {
                bits = (int) ((binaryKmer[arrayIndex] >>> (62 - bitOffset)) & 0b11);
            } else {
                int bitsInFirstLong = 64 - bitOffset;
                long part1 = (binaryKmer[arrayIndex] & ((1L << bitsInFirstLong) - 1)) << (2 - bitsInFirstLong);
                long part2 = (binaryKmer[arrayIndex + 1] >>> (64 - (2 - bitsInFirstLong)));
                bits = (int) ((part1 | part2) & 0b11);
            }

            // get complement
            int compBits = (~bits) & 0b11;

            // set in reverse array
            int revIndex = (kmerLength - i - 1) * 2;
            int revArrayIndex = revIndex / 64;
            int revOffset = revIndex % 64;

            if (revOffset <= 62) {
                reverse[revArrayIndex] |= ((long) compBits) << (62 - revOffset);
            } else {
                int bitsInFirstLong = 64 - revOffset;
                reverse[revArrayIndex] |= ((long) compBits) >>> (2 - bitsInFirstLong);
                reverse[revArrayIndex + 1] |= ((long) compBits) << (62 - (2 - bitsInFirstLong));
            }
        }

        return reverse;
    }

    /***
     * convert the kmer long array to a string. used in the toString method to convert the kmer long array to a string.
     * @param binaryKmer
     * @param kmerLength
     * @return
     */
    public static String longToKmer(long binaryKmer, int kmerLength) {
        return longToKmer(new long[]{binaryKmer}, kmerLength);
    }

    /***
     * convert the kmer long array to a string. used in the toString method to convert the kmer long array to a string.
     * @param binaryKmer
     * @param kmerLength
     * @return
     */
    public static String longToKmer(long[] binaryKmer, int kmerLength) {
        char[] kmer = new char[kmerLength];

        for (int i = 0; i < kmerLength; i++) {
            int bitIndex = i * 2;
            int arrayIndex = bitIndex / 64;
            int bitOffset = bitIndex % 64;

            int bits;
            if (bitOffset <= 62) {
                bits = (int) ((binaryKmer[arrayIndex] >>> (62 - bitOffset)) & 0b11);
            } else {
                int bitsInFirstLong = 64 - bitOffset;
                long part1 = (binaryKmer[arrayIndex] & ((1L << bitsInFirstLong) - 1)) << (2 - bitsInFirstLong);
                long part2 = (binaryKmer[arrayIndex + 1] >>> (64 - (2 - bitsInFirstLong)));
                bits = (int) ((part1 | part2) & 0b11);
            }

            kmer[i] = bitsToBase(bits);
        }

        return new String(kmer);
    }

    /***
     * Converts a binary representation of a k-mer to a nucleotide string.
     */
    public static char[] binaryToKmer(long binaryKmer, int kmerLength) {
        char[] kmer = new char[kmerLength];
        for (int i = kmerLength - 1; i >= 0; i--) {
            kmer[i] = bitsToBase((int) (binaryKmer & 0b11)); // Map the last two bits to a base
            binaryKmer >>= 2; // Shift right to process the next base
        }
        return kmer;
    }

    /***
     * convert the bits to base. used in the longToKmer method to convert the kmer long array to a string.
     * @param twoBits
     * @return
     */
    private static char bitsToBase(int twoBits) {
        return switch (twoBits) {
            case 0b00 -> 'A';
            case 0b01 -> 'C';
            case 0b10 -> 'G';
            case 0b11 -> 'T';
            default -> throw new IllegalArgumentException("Invalid bits for base: " + twoBits);
        };
    }

    /***
     * compare two long arrays. used in the compareTo method to compare two kmer long arrays.
     * @param a
     * @param b
     * @return
     */
    private static int compareLongArrays(long[] a, long[] b) {
        int len = Math.min(a.length, b.length);
        for (int i = 0; i < len; i++) {
            if (a[i] != b[i]) {
                return Long.compareUnsigned(a[i], b[i]);
            }
        }
        return Integer.compare(a.length, b.length);
    }

    @Override
    public int compareTo(Kmer o) {
        return compareLongArrays(this.kmerLong, o.kmerLong);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Kmer) {
            Kmer other = (Kmer) obj;
            return Arrays.equals(this.kmerLong, other.kmerLong);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(this.kmerLong);
    }

    @Override
    public String toString() {
        return longToKmer(this.kmerLong, this.kmerLength);
    }

    /***
     * get the prefix forward as a string. used in the toString method to convert the kmer long array to a string. [miscellaneous]
     * @return
     */
    public String getPrefixFwdStr() {
        return longToKmer(getPrefixFwd(), this.prefixLength);
    }

    /***
     * TODO: this could improve performance
     * get the suffix forward as a string. used in the toString method to convert the kmer long array to a string. [miscellaneous]
     * @param base
     * @return
     */
    public Kmer insertBase(char base) {
        if (this.isCanonical) {
            Logger.error(CLASS_NAME, "Cannot insert base into canonical kmer");
            throw new IllegalArgumentException("Cannot insert base into canonical kmer");
        }

        long[] newKmerLong = new long[this.kmerLong.length];

        // shift left and insert base bits
        for (int i = 0; i < this.kmerLong.length; i++) {
            if (this.kmerLong.length == (i + 1)) {
                newKmerLong[i] = (this.kmerLong[i] << 2) | baseToBits(base);
            } else {
                newKmerLong[i] = (this.kmerLong[i] << 2) | ((this.kmerLong[i + 1] >> 62) & 0b11);
            }
        }

        // mask unused bits in the most significant word
        int totalBits = 2 * this.kmerLength;
        int unusedBits = (64 * this.kmerLong.length) - totalBits;
        if (unusedBits > 0) {
            long mask = -1L >>> unusedBits;  // Keep only the low `totalBits` bits
            newKmerLong[0] &= mask;          // Apply mask to the highest word
        }

        return new Kmer(newKmerLong, this.kmerLength, this.prefixLength, this.isCanonical);
    }

}
//EOF