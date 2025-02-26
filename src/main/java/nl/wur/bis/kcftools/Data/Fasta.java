package nl.wur.bis.kcftools.Data;

import java.io.BufferedWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import static nl.wur.bis.kcftools.Utils.HelperFunctions.fold_seq;

/***
 * This class represents a Fasta object with thread-safe operations.
 */
public class Fasta implements Serializable {
    private static final int FOLDLEN = 60;
    private final Object lock = new Object(); // Synchronization lock for thread safety

    private int id;
    private String name;
    private String sequence;
    private String description;

    // no-arg constructor
    public Fasta() {
    }

    public Fasta(int id, String name, String sequence) {
        this.id = id;
        this.name = name;
        this.sequence = sequence;
        this.description = "";
    }

    public Fasta(int id, String name, String sequence, String description) {
        this.id = id;
        this.name = name;
        this.sequence = sequence;
        this.description = description;
    }

    public int getId() {
        synchronized (lock) {
            return id;
        }
    }

    public String getName() {
        synchronized (lock) {
            return name;
        }
    }

    public String getSequence() {
        synchronized (lock) {
            return sequence;
        }
    }

    public String getQuality() {
        return ""; // No mutable state, remains safe
    }

    public String getDescription() {
        synchronized (lock) {
            return description;
        }
    }

    public void setName(String name) {
        synchronized (lock) {
            this.name = name;
        }
    }

    public void setSequence(String sequence) {
        synchronized (lock) {
            this.sequence = sequence;
        }
    }

    public void setQuality(String s) {
        // do nothing
    }

    public void setId(int id) {
        synchronized (lock) {
            this.id = id;
        }
    }

    public void writeFasta(BufferedWriter writer) {
        synchronized (lock) {
            try {
                writer.write(">" + name + " " + description + "\n");
                writer.write(fold_seq(sequence, FOLDLEN));
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void setDescription(String description) {
        synchronized (lock) {
            this.description = description;
        }
    }

    @Override
    public String toString() {
        synchronized (lock) {
            return ">" + name + " " + description + "\n" + fold_seq(sequence, FOLDLEN);
        }
    }

    public String toFastaString(boolean fold) {
        synchronized (lock) {
            if (fold) {
                return ">" + name + " " + description + "\n" + fold_seq(sequence, FOLDLEN);
            } else {
                return ">" + name + " " + description + "\n" + sequence + "\n";
            }
        }
    }

    public int length() {
        synchronized (lock) {
            return sequence.length();
        }
    }

    public int getLength() {
        synchronized (lock) {
            return sequence.length();
        }
    }

    public String getUppercaseSequence() {
        synchronized (lock) {
            return sequence.toUpperCase();
        }
    }

    /***
     * This method will return the kmers of the sequence as a list (NOT UNIQUE KMERS)
     */
    public List<Kmer> getKmersList(int kmerLength, int prefixLength, boolean useCanonical) {
        synchronized (lock) {
            List<Kmer> mers = new ArrayList<>();
            char[] kmerChars = new char[kmerLength];
            Kmer kmer = null;

            int validStart = 0;
            for (int i = 0; i < sequence.length(); i++) {
                char base = Character.toUpperCase(sequence.charAt(i));
                if (!isValidBase(base)) {
                    validStart = i + 1;
                    kmer = null;
                    continue;
                }

                int offset = i - validStart;
                if (offset < kmerLength) {
                    kmerChars[offset] = base;
                    if (offset == kmerLength - 1) {
                        kmer = new Kmer(kmerChars, prefixLength, useCanonical);
                    }
                } else if (kmer != null) {
                    kmer = kmer.insertBase(base);
                }

                if (kmer != null) {
                    mers.add(kmer);
                }
            }
            return mers;
        }
    }

    /***
     * This method will return the true if the base is valid
     */
    private boolean isValidBase(char base) {
        return base == 'A' || base == 'C' || base == 'G' || base == 'T';
    }

    public int getEffectiveATGCCount(int kmerLength) {
        synchronized (lock) {
            int count = 0;
            int stretchLength = 0; // Track the length of consecutive ATGC stretch

            for (int i = 0; i < sequence.length(); i++) {
                char base = Character.toUpperCase(sequence.charAt(i));

                // Check if the base is one of A, T, G, C
                if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                    stretchLength++;
                } else {
                    // End of a stretch: check if it's long enough
                    if (stretchLength >= kmerLength) {
                        count += stretchLength;
                    }
                    stretchLength = 0; // Reset stretch length
                }
            }

            // Handle the last stretch, if the sequence ends with ATGC
            if (stretchLength >= kmerLength) {
                count += stretchLength;
            }

            return count;
        }
    }
}
