package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.Objects;

public class Kmer implements Comparable<Kmer> {
    private final int kmerLength;
    private final int prefixLength;
    private final int suffixLength;
    private long kmerLong = 0;
    private final boolean isCanonical;

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
        long kmerLong;
        this.kmerLength = kmer.length;
        this.prefixLength = prefixLength;
        this.suffixLength = this.kmerLength - prefixLength;
        this.isCanonical = isCanonical;
        this.kmerLong = kmerToLong(kmer);
        if (isCanonical) {
            // get the lowest of the kmer and its reverse complement
            this.kmerLong = Math.min(this.kmerLong, getReverseComplement(this.kmerLong, this.kmerLength));
        }
    }


    public Kmer(long kmerLong, int kmerLength, int prefixLength, boolean isCanonical) {
        this.kmerLength = kmerLength;
        this.prefixLength = prefixLength;
        this.suffixLength = this.kmerLength - prefixLength;
        this.isCanonical = isCanonical;
        this.kmerLong = kmerLong;
        if (isCanonical) {
            // get the lowest of the kmer and its reverse complement
            this.kmerLong = Math.min(this.kmerLong, getReverseComplement(this.kmerLong, this.kmerLength));
        }
    }

    public Kmer(Kmer kmer, boolean isCanonical){
        // copy the kmer and set the canonical status
        this.kmerLength = kmer.getKmerLength();
        this.prefixLength = kmer.getPrefixLength();
        this.suffixLength = kmer.getSuffixLength();
        this.isCanonical = isCanonical;
        this.kmerLong = kmer.getKmerLong();
        if (isCanonical) {
            // get the lowest of the kmer and its reverse complement
            this.kmerLong = Math.min(this.kmerLong, getReverseComplement(this.kmerLong, this.kmerLength));
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

    public long getKmerLong() {
        return this.kmerLong;
    }

    /***
     * Lazily calculates and returns the forward prefix as an integer.
     */
    public int getPrefixFwd() {
        if (this.prefixFwd == null) {
            this.prefixFwd = (int) (this.kmerLong >> this.suffixLength * 2);
        }
        return this.prefixFwd;
    }

    /***
     * Lazily calculates and returns the forward suffix as a byte array.
     */
    public byte[] getSuffixFwd() {
        if (this.suffixFwd == null) {
            this.suffixFwd = extractSuffix(this.kmerLong, this.suffixLength);
        }
        return this.suffixFwd;
    }


    /***
     * Lazily calculates and returns the reverse prefix as an integer.
     */
    public int getPrefixRev() {
        if (this.prefixRev == null) {
            // reverse complement the kmer and get the prefix
            this.prefixRev = (int) (getReverseComplement(this.kmerLong, this.kmerLength) >> this.suffixLength * 2);
        }
        return this.prefixRev;
    }

    /***
     * Lazily calculates and returns the reverse suffix as a byte array.
     */
    public byte[] getSuffixRev() {
        if (this.suffixRev == null) {
            // reverse complement the kmer and get the suffix
            this.suffixRev = extractSuffix(getReverseComplement(this.kmerLong, this.kmerLength), this.suffixLength);
        }
        return this.suffixRev;
    }

    /***
     * Converts the kmer's binary representation back to a nucleotide string.
     * @return Nucleotide string
     */
    @Override
    public String toString() {
        return longToKmer(this.kmerLong, this.kmerLength);
    }

    @Override
    public int compareTo(Kmer o) {
        return Long.compare(this.kmerLong, o.kmerLong);
    }

    /***
     * Check if kmer is equal to another kmer.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Kmer) {
            Kmer other = (Kmer) obj;
            return this.kmerLong == other.kmerLong;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(this.kmerLong);
    }

    /***
     * Generates the forward suffix as a nucleotide string.
     * @return Suffix nucleotide string
     */
    public String getSuffixForwardString() {
        return convertBytesToKmer(getSuffixFwd(), this.suffixLength);
    }

    /***
     * Generates the forward prefix as a nucleotide string.
     * @return Prefix nucleotide string
     */
    public String getPrefixForwardString() {
        return longToKmer(getPrefixFwd(), this.prefixLength);
    }

    /***
     * Helper function to extract suffix from a binary representation.
     * @param binaryKmer The binary representation of the k-mer
     * @param suffixLength The length of the suffix to extract
     * @return Byte array representation of the suffix
     */
    private byte[] extractSuffix(long binaryKmer, int suffixLength) {
        byte[] suffix = new byte[(suffixLength + 3) / 4];
        long suffixBits = binaryKmer & ((1L << (2 * suffixLength)) - 1);

        // Pack the suffix into the byte array
        for (int i = 0; i < suffixLength; ++i) {
            int shift = (suffixLength - i - 1) * 2;
            suffix[i / 4] |= (byte) (((suffixBits >> shift) & 0x3) << ((3 - i % 4) * 2));
        }
        return suffix;
    }

    /***
     * Converts a byte array representation of a suffix to its nucleotide string.
     * @param suffixBytes The byte array containing suffix information
     * @param suffixLength The length of the suffix
     * @return Suffix nucleotide string
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
            default -> {
                Logger.error(CLASS_NAME, "Invalid base index: " + baseIndex);
                throw new IllegalArgumentException("Invalid base index: " + baseIndex);
            }
        };
    }

    /***
     * This method is used to get the signature of a kmer
     */
    public int getSignature(Signature signature) {

        // get length of the signature to be computed form signature object
        int signatureLength = signature.getSignLength();
        // get the first signature length bases of the kmer
        int currentSignature = (int) (this.kmerLong >> (2 * (this.kmerLength - signatureLength)));
        // Constrain currentSignature to valid range
        int minSignature = signature.getSignature(currentSignature);
        for (int i = signatureLength; i < this.kmerLength; i++) {
            currentSignature = (int) ((currentSignature << 2) | ((this.kmerLong >> (2 * (this.kmerLength - i - 1))) & 0x03));
            // Constrain currentSignature to valid range
            currentSignature &= (1 << (2 * signatureLength)) - 1;


            // Check if the current signature is within the allowed range
            if (currentSignature >= (1 << (2 * signatureLength))) {
                Logger.error(CLASS_NAME, "Invalid signature (out of range): " + currentSignature);
            }

            // Check if the signature is valid
            if (signature.getSignature(currentSignature) < minSignature) {
                minSignature = signature.getSignature(currentSignature);
            }
        }

        return minSignature;
    }

    public static long kmerToLong(String kmer) {
        return kmerToLong(kmer.toCharArray());
    }

    /***
     * This method is used to get the canonical status of a kmer
     */
    public static long kmerToLong(char[] kmer) {
        long binaryKmer = 0;
        for (char base : kmer) {
            binaryKmer <<= 2; // Shift left by 2 bits
            binaryKmer |= baseToBits(base); // Add the bits for the current base
        }
        return binaryKmer;
    }

    /***
     * Converts a nucleotide base to a two-bit representation.
     */
    private static int baseToBits(char base) {
        return switch (base) {
            case 'A' -> 0b00; // A -> 00
            case 'C' -> 0b01; // C -> 01
            case 'G' -> 0b10; // G -> 10
            case 'T' -> 0b11; // T -> 11
            default -> throw new IllegalArgumentException("Invalid nucleotide in k-mer: " + base);
        };
    }

    /***
     * Converts a binary representation of a k-mer to a nucleotide string.
     */
    public static String longToKmer(long binaryKmer, int kmerLength) {
        return new String(binaryToKmer(binaryKmer, kmerLength));
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
     * Converts a two-bit representation of a base to its character representation.
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
     * This method is used to get the reverse complement of a kmer in binary form
     */
    public static long getReverseComplement(long binaryKmer, int kmerLength) {
        long reverseComplement = 0;
        for (int i = 0; i < kmerLength; i++) {
            // Extract the last 2 bits (base)
            long base = binaryKmer & 0b11;
            // Complement the base: A (00) -> T (11), C (01) -> G (10), G (10) -> C (01), T (11) -> A (00)
            base = (~base) & 0b11;
            // Shift the complemented base to the reverse position
            reverseComplement = (reverseComplement << 2) | base;
            // Shift the original binaryKmer to process the next base
            binaryKmer >>= 2;
        }
        return reverseComplement;
    }


    public String getPrefixReverseString() {
        return longToKmer(getPrefixRev(), this.prefixLength);
    }

    public String getSuffixReverseString() {
        return convertBytesToKmer(getSuffixRev(), this.suffixLength);
    }

    public String getRevKmerString() {
        if (this.reverseComplementKmer == null) {
            this.reverseComplementKmer = longToKmer(getReverseComplement(this.kmerLong, this.kmerLength), this.kmerLength);
        }
        return this.reverseComplementKmer;
    }

    /***
     * this method will insert a base at the end of the kmer and remove the first base and construct a new kmer
     * make sure to use this with caution of the kmer is not canonical
     */
    public Kmer insertBase(char base) {
        if (this.isCanonical) {
            Logger.error(CLASS_NAME, "Cannot insert base into canonical kmer");
            throw new IllegalArgumentException("Cannot insert base into canonical kmer");
        }

        // Convert the base to its 2-bit representation
        long newBaseBits = baseToBits(base);

        // Left-shift kmerLong by 2 bits to make space for the new base
        long newKmer = (this.kmerLong << 2) | newBaseBits;

        // Mask to ensure newKmer remains within the bounds of the k-mer length
        long mask = (1L << (2 * this.kmerLength)) - 1; // Create a mask for the k-mer length
        newKmer &= mask;

        return new Kmer(newKmer, this.kmerLength, this.prefixLength, this.isCanonical);
    }
}
//EOF