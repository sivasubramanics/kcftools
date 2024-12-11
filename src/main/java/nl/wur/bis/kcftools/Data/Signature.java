package nl.wur.bis.kcftools.Data;


/***
 * This is a class to initialise and process the signatures of the kmers
 * inspired by the KMC2 implementation
 */
public class Signature {

    private final int signLength;
    private final int[] signRefMap;

    public Signature(int signLength) {
        this.signLength = signLength;
        this.signRefMap = new int[1 << (signLength * 2)];
        initNorm();
    }

    public int getSignature(int mmer) {
        return signRefMap[mmer];
    }

    private void initNorm() {
        int special = 1 << (signLength * 2); // Calculate the special value
        for (int i = 0; i < special; ++i) {
            // Get the reverse complement of the current value
            int rev = getRev(i, signLength);

            // Determine the value for the forward and reverse strand
            int strVal = isAllowed(i, signLength) ? i : special;
            int revVal = isAllowed(rev, signLength) ? rev : special;

            // Store the minimum of the two in the norm array
            signRefMap[i] = Math.min(strVal, revVal);
        }
    }

    /***
     * Check if the signature is allowed (inspired by the KMC2 implementation)
     */
    private static boolean isAllowed(int signature, int signLength) {
        // Check for disallowed suffixes
        if ((signature & 0x3F) == 0x3F) { // TTT suffix
            return false;
        }
        if ((signature & 0x3F) == 0x3B) { // TGT suffix
            return false;
        }
        if ((signature & 0x3C) == 0x3C) { // TG* suffix
            return false;
        }

        // Check for disallowed patterns inside the signature
        for (int j = 0; j < signLength - 3; ++j) {
            if ((signature & 0xF) == 0) { // AA inside
                return false;
            } else {
                signature >>= 2; // Shift right to check the next 2 bits
            }
        }

        // Check for disallowed prefixes
        if (signature == 0) { // AAA prefix
            return false;
        }
        if (signature == 0x04) { // ACA prefix
            return false;
        }
        if ((signature & 0xF) == 0) { // *AA prefix
            return false;
        }

        return true; // If none of the disallowed patterns are found, return true
    }


    /***
     * Get the reverse complement of a sequence
     * @param sequence The sequence to get the reverse complement of
     * @param length The length of the sequence
     * @return The reverse complement of the sequence
     */
    private int getRev(int sequence, int length) {
        int reverseComplement = 0;
        for (int i = 0; i < length; i++) {
            // Get the last 2 bits (base)
            int base = sequence & 0b11;
            // Complement the base (A=0->T=3, C=1->G=2, G=2->C=1, T=3->A=0)
            base = (~base) & 0b11;
            // Shift the complemented base to the reverse position
            reverseComplement = (reverseComplement << 2) | base;
            // Shift the sequence right to process the next base
            sequence >>= 2;
        }
        return reverseComplement;
    }

    public int getSignLength() {
        return signLength;
    }
}
//EOF