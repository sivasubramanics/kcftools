package nl.wur.bis.kcftools.Data;

import java.io.*;
import java.util.*;
import nl.wur.bis.kcftools.Utils.Logger;
import org.jgrapht.*;
import org.jgrapht.graph.*;

/**
 * Class to parse GTF files.
 * It builds a directed graph of features and allows exporting to GTF format.
 * Works well and tested with AGAT formatted GTF files.
 * Highly recommended to use AGAT to format GTF files.
 */
public class GTF {
    private String filePath;
    private Graph<String, DefaultEdge> features;
    private HashMap<String, Feature> featureMap;
    private final String CLASS_NAME = this.getClass().getSimpleName();

    public GTF(String filePath) {
        this.filePath = filePath;
        parseGTF();
    }

    private void parseGTF() {
        Logger.info(CLASS_NAME, "Parsing GTF file at: " + filePath);
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            HashMap<String, Integer> exonCounts = new HashMap<>();
            features = new DefaultDirectedGraph<>(DefaultEdge.class);
            featureMap = new HashMap<>();

            while ((line = br.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;
                String[] fields = line.split("\t");
                if (fields.length < 9) {
                    Logger.error(CLASS_NAME, "Malformed line: " + line);
                    continue;
                }

                HashMap<String, String> attributes = parseAttributes(fields[8]);
                String type = fields[2];
                String featureID = null, parentID = null, chromID = fields[0];

                features.addVertex(chromID); // always add chromosome node

                switch (type) {
                    case "gene": case "pseudogene":
                        featureID = attributes.get("gene_id");
                        parentID = chromID;
                        break;
                    case "transcript": case "mRNA": case "RNA": case "lnc_RNA":
                    case "rRNA": case "tRNA": case "snRNA": case "snoRNA":
                        featureID = attributes.get("transcript_id");
                        parentID = attributes.get("gene_id");
                        if (featureID.equals(parentID)) {
                            Logger.error(CLASS_NAME, "Transcript ID is the same as Gene ID: " + featureID + ". Fix the GTF file using AGAT.");
                            continue;
                        }
                        // in case the GTF file does not have gene feature, we create a gene feature
                        if (!features.containsVertex(parentID)) {
                            features.addVertex(parentID);
                            Feature geneFeature = new Feature(chromID, Integer.parseInt(fields[3]), Integer.parseInt(fields[4]),
                                    fields[6].charAt(0), "gene", parentID);
                            features.addEdge(chromID, parentID);
                            featureMap.put(parentID, geneFeature);
                        }
                        // update start and end of the gene feature in featureMap and assign them back to the featureMap
                        Feature geneFeature = featureMap.get(parentID);
                        if (geneFeature != null) {
                            geneFeature.updateStart(Integer.parseInt(fields[3]));
                            geneFeature.updateEnd(Integer.parseInt(fields[4]));
                            featureMap.put(parentID, geneFeature);
                        }
                        break;
                    case "exon":
                        parentID = attributes.get("transcript_id");
                        int count = exonCounts.getOrDefault(parentID, 0) + 1;
                        exonCounts.put(parentID, count);
                        featureID = parentID + "-e-" + count;
                        break;
                    default:
                        continue;
                }

                Feature feature = new Feature(fields[0], Integer.parseInt(fields[3]), Integer.parseInt(fields[4]),
                        fields[6].charAt(0), type, featureID);
                features.addVertex(featureID);
                featureMap.put(featureID, feature);

                if (parentID != null) {
                    features.addVertex(parentID);
                    features.addEdge(parentID, featureID);
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error parsing GTF file: " + e.getMessage());
        }
    }

    /**
     * TEST UTILITY FUNCTION
     * Exports the GTF features to a file in GTF format.
     * The output will contain gene, transcript, and exon features.
     * The source is set to "KCFtools".
     */
    public void exportGTF(String outputFile) {
        final String SOURCE = "KCFtools";

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            for (String chrom : getChromosomes()) {
                for (String geneID : getGenes(chrom)) {
                    writeFeature(writer, featureMap.get(geneID), featureMap.get(geneID).type, geneID, null, null, SOURCE);

                    for (String transcriptID : getTranscripts(geneID)) {
                        writeFeature(writer, featureMap.get(transcriptID), featureMap.get(transcriptID).type,
                                geneID, transcriptID, null, SOURCE);

                        for (String exonID : getExons(transcriptID)) {
                            writeFeature(writer, featureMap.get(exonID), "exon",
                                    geneID, transcriptID, exonID, SOURCE);
                        }
                    }
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Failed to export GTF to file: " + outputFile + ". Reason: " + e.getMessage());
        }
    }

    /**
     * TEST UTILITY FUNCTION
     * Writes a feature to the GTF file.
     * The feature can be a gene, transcript, or exon.
     */
    private void writeFeature(BufferedWriter writer, Feature feature, String featureType,
                              String geneID, String transcriptID, String exonID, String source) throws IOException {
        if (feature == null) return;

        StringBuilder attributes = new StringBuilder();
        attributes.append(String.format("gene_id \"%s\";", geneID));
        if (transcriptID != null) attributes.append(String.format(" transcript_id \"%s\";", transcriptID));
        if (exonID != null) attributes.append(String.format(" exon_id \"%s\";", exonID));

        String line = String.format("%s\t%s\t%s\t%d\t%d\t.\t%c\t.\t%s\n",
                feature.chromosome, source, featureType, feature.start, feature.end, feature.strand, attributes.toString());

        writer.write(line);
    }

    /**
     * Parses the attributes string from a GTF line.
     * The attributes are expected to be in the format: key "value"; key "value"; ...
     */
    private HashMap<String, String> parseAttributes(String attrStr) {
        HashMap<String, String> map = new HashMap<>();
        for (String attr : attrStr.split(";")) {
            String[] pair = attr.trim().replaceAll("\"", "").split(" ");
            if (pair.length == 2) map.put(pair[0], pair[1]);
        }
        return map;
    }

    /**
     * Returns an array of chromosome IDs present in the GTF file.
     * Chromosomes are identified as vertices with no incoming edges.
     */
    public String[] getChromosomes() {
        List<String> chroms = new ArrayList<>();
        for (String v : features.vertexSet()) {
            if (features.inDegreeOf(v) == 0) {
                chroms.add(v);
            }
        }
        return chroms.toArray(new String[0]);
    }

    /**
     * Returns an array of gene IDs for a given chromosome ID.
     * Genes are identified as vertices with chromosome as parent
     */
    public String[] getGenes(String chromID) {
        return getChildren(chromID);
    }

    /**
     * Returns an array of transcript IDs for a given gene ID.
     * Transcripts are identified as vertices with gene as parent.
     */
    public String[] getTranscripts(String geneID) {
        return getChildren(geneID);
    }

    /**
     * Returns an array of exon IDs for a given transcript ID.
     * Exons are identified as vertices with transcript as parent.
     */
    public String[] getExons(String transcriptID) {
        return getChildren(transcriptID);
    }

    /**
     * Returns an array of children IDs for a given parent ID.
     * Children are identified as vertices with parent as source.
     */
    private String[] getChildren(String parentID) {
        if (!features.containsVertex(parentID)) return new String[0];
        List<String> children = new ArrayList<>();
        for (DefaultEdge e : features.outgoingEdgesOf(parentID)) {
            String child = features.getEdgeTarget(e);
            if (!child.equals(parentID)) {
                children.add(child);
            }
        }
        return children.toArray(new String[0]);
    }

    /**
     * Returns a Fasta object for a given feature ID.
     * The Fasta object contains the sequence of the feature and its description.
     */
    public Fasta getFasta(String featureID, FastaIndex fastaIndex, boolean isGene) {
        if (!features.containsVertex(featureID)) return null;
        Set<GTF.Loci> lociSet = new HashSet<>();
        String[] targets = isGene ? getTranscripts(featureID) : getExons(featureID);

        for (String t : targets) {
            String[] exons = isGene ? getExons(t) : new String[]{t};
            for (String exonID : exons) {
                Feature f = featureMap.get(exonID);
                if (f != null)
                    lociSet.add(new Loci(f.chromosome, f.start, f.end, String.valueOf(f.strand)));
            }
        }
        if (lociSet.isEmpty()) return null;

        List<Loci> merged = mergeOverlappingLoci(lociSet);
        Collections.sort(merged);
        StringBuilder seq = new StringBuilder();
        for (Loci loci : merged) {
            String s = fastaIndex.getSequence(loci.chromosome, loci.start - 1, loci.getLength());
            if (s != null) seq.append(s);
        }

        String desc = String.join(" ", merged.stream().map(l -> l.chromosome + ":" + l.start + "-" + l.end + "[" + l.strand + "]").toList());
        return new Fasta(-1, featureID, seq.toString(), desc);
    }

    /**
     * Writes the FASTA sequences for genes or transcripts to a file.
     * If isGene is true, it writes gene sequences; otherwise, it writes transcript sequences.
     */
    public void writeFasta(String outputFile, FastaIndex fastaIndex, boolean isGene) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            for (String chrom : getChromosomes()) {
                for (String geneID : getGenes(chrom)) {
                    if (isGene) {
                        Fasta f = getFasta(geneID, fastaIndex, true);
                        if (f != null) f.writeFasta(writer);
                    } else {
                        for (String tID : getTranscripts(geneID)) {
                            Fasta f = getFasta(tID, fastaIndex, false);
                            if (f != null) f.writeFasta(writer);
                        }
                    }
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error writing FASTA file: " + e.getMessage());
        }
    }

    /***
     * Merges overlapping loci from a set of Loci objects.
     * Overlapping loci are merged into a single Loci object.
     */
    public static List<Loci> mergeOverlappingLoci(Set<Loci> lociSet) {
        List<Loci> sorted = new ArrayList<>(lociSet);
        Collections.sort(sorted);
        List<Loci> merged = new ArrayList<>();
        for (Loci current : sorted) {
            if (merged.isEmpty()) merged.add(current);
            else {
                Loci last = merged.get(merged.size() - 1);
                if (last.overlapsWith(current))
                    merged.set(merged.size() - 1, last.mergeWith(current));
                else
                    merged.add(current);
            }
        }
        return merged;
    }

    /***
     * Returns the Loci object for a given feature ID.
     * The Loci object contains the chromosome, start, end, and strand of the feature.
     */
    public Loci getLoci(String featureID) {
        if (!featureMap.containsKey(featureID)) {
            Logger.error(CLASS_NAME, "Feature ID not found: " + featureID);
            return null;
        }
        Feature feature = featureMap.get(featureID);
        return new Loci(feature.chromosome, feature.start, feature.end, String.valueOf(feature.strand));
    }

    /***
     * SubClass representing a feature in the GTF file.
     */
    public static class Feature implements Comparable<Feature> {
        private final String chromosome;
        private int start;
        private int end;
        private char strand;
        private final String type;
        private final String id;

        public Feature(String chromosome, int start, int end, char strand, String type, String id) {
            this.chromosome = chromosome;
            this.start = start;
            this.end = end;
            this.strand = strand;
            this.type = type;
            this.id = id;
        }

        @Override
        public int compareTo(Feature other) {
            if (this.chromosome.equals(other.chromosome)) {
                return Integer.compare(this.start, other.start);
            }
            return this.chromosome.compareTo(other.chromosome);
        }

        @Override
        public String toString() {
            return id + "\t" + type + "\t" + chromosome + "\t" + start + "\t" + end + "\t" + strand;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Feature)) return false;
            Feature other = (Feature) o;
            return this.id.equals(other.id);
        }

        @Override
        public int hashCode() {
            return id.hashCode();
        }

        public String getLoci() {
            return new Loci(chromosome, start, end, String.valueOf(strand)).toString();
        }

        public void updateStart(int start) {
            if (start < this.start) this.start = start;
        }

        public void updateEnd(int end) {
            if (end > this.end) this.end = end;
        }

    }

    /***
     * SubClass representing a genomic loci with chromosome, start, end, and strand.
     * Implements Comparable for sorting and merging functionalities.
     */
    public static class Loci implements Comparable<Loci> {
        private final String chromosome;
        private final int start;
        private final int end;
        private final String strand;

        public Loci(String chromosome, int start, int end, String strand) {
            this.chromosome = chromosome;
            this.start = start;
            this.end = end;
            this.strand = strand;
        }

        public String getChromosome() {
            return chromosome;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public String getStrand() {
            return strand;
        }

        public boolean overlapsWith(Loci other) {
            if (!this.chromosome.equals(other.chromosome)) return false;
            if (!this.strand.equals(other.strand)) return false;
            return this.start <= other.end && other.start <= this.end;
        }

        public Loci mergeWith(Loci other) {
            int mergedStart = Math.min(this.start, other.start);
            int mergedEnd = Math.max(this.end, other.end);
            return new Loci(this.chromosome, mergedStart, mergedEnd, this.strand);
        }

        @Override
        public int compareTo(Loci other) {
            if (this.chromosome.equals(other.chromosome)) {
                return Integer.compare(this.start, other.start);
            }
            return this.chromosome.compareTo(other.chromosome);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Loci)) return false;
            Loci other = (Loci) o;
            return this.chromosome.equals(other.chromosome)
                    && this.start == other.start
                    && this.end == other.end
                    && this.strand.equals(other.strand);
        }

        @Override
        public int hashCode() {
            int result = chromosome.hashCode();
            result = 31 * result + Integer.hashCode(start);
            result = 31 * result + Integer.hashCode(end);
            result = 31 * result + strand.hashCode();
            return result;
        }

        public int getLength() {
            return end - start + 1;
        }
    }
}
// EOF