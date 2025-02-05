package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.HelperFunctions;

import java.io.*;
import java.util.*;
import java.util.stream.Stream;

import static nl.wur.bis.kcftools.Utils.HelperFunctions.reverseComplement;

public class GTFReader {

    private Map<String, Feature> genes;
    private Map<String, String> transcriptToGene;
    private ArrayList<String> geneIDs;
    private ArrayList<String> transcriptIDs;

    public GTFReader() {
    }

    public GTFReader(String filePath) throws IOException {
        genes = new LinkedHashMap<>();
        geneIDs = new ArrayList<>();
        transcriptIDs = new ArrayList<>();
        transcriptToGene = new HashMap<>();
        readGTF(filePath);
    }

    // Reads a GTF file and builds the hierarchical structure
    public void readGTF(String filePath) throws IOException {
        HelperFunctions.log("info", "GTFReader", "Reading GTF file: " + filePath);
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue; // Skip comments

                String[] fields = line.split("\t");
                if (fields.length < 9) continue;

                String chrom = fields[0];
                String type = fields[2];
                int start = Integer.parseInt(fields[3]);
                int end = Integer.parseInt(fields[4]);
                String strand = fields[6];
                String attributes = fields[8];

                // Extract identifiers
                String geneId = extractAttribute(attributes, "gene_id");
                String transcriptId = extractAttribute(attributes, "transcript_id");

                // Create a Feature object
                String featureId = type.equals("gene") ? geneId : transcriptId;
                if (featureId == null) continue;

                Feature feature = new Feature(type, chrom, start, end, strand, featureId);

                // Link features based on type
                if (type.equals("gene")) {
                    genes.putIfAbsent(geneId, feature);
                    geneIDs.add(geneId);
                } else if (type.equals("transcript")) {
                    Feature gene = genes.get(geneId);
                    if (gene != null) {
                        feature.parent = gene;
                        gene.addChild(feature);
                        transcriptToGene.put(transcriptId, geneId); // Map transcript to its parent gene
                        transcriptIDs.add(transcriptId);
                    }
                } else { // For exons, CDS, etc.
                    String parentGeneId = transcriptToGene.get(transcriptId);
                    if (parentGeneId != null) {
                        Feature gene = genes.get(parentGeneId);
                        if (gene != null) {
                            for (Feature transcript : gene.children) {
                                if (transcript.id.equals(transcriptId)) {
                                    feature.parent = transcript;
                                    transcript.addChild(feature);
                                }
                            }
                        }
                    }
                }
            }

            // Sort all children by start position
            for (Feature gene : genes.values()) {
                gene.sortChildrenByStart();
                for (Feature transcript : gene.children) {
                    transcript.sortChildrenByStart();
                }
            }
        }

        HelperFunctions.log("info", "GTFReader", "Finished reading GTF file: " + filePath);
        if (genes.isEmpty()) {
            HelperFunctions.log("warn", "GTFReader", "No 'gene' features found in GTF file: " + filePath);
            HelperFunctions.log("warn", "GTFReader", "Please check the GTF file format, else fix them using agat");
            HelperFunctions.log("error", "GTFReader", "Exiting the program");
        }
    }

    /**
     * Extract an attribute value from the attributes string
     */
    private String extractAttribute(String attributes, String key) {
        String[] pairs = attributes.split(";");
        for (String pair : pairs) {
            String[] kv = pair.trim().split(" ");
            if (kv.length == 2 && kv[0].equals(key)) {
                return kv[1].replaceAll("\"", ""); // Remove quotes
            }
        }
        return null;
    }

    /**
     * Get fasta sequence for a feature by ID and type
     */
    public Fasta getFasta(String id, String featureType, FastaIndex fastaIndex) {
        return switch (featureType.toLowerCase()) {
            case "gene" -> getGeneFasta(id, fastaIndex);
            case "transcript" -> getTranscriptFasta(id, fastaIndex);
            case "exon" -> getExonFasta(id, fastaIndex);
            case "cds" -> getCDSFasta(id, fastaIndex);
            default -> null;
        };
    }

    /**
     * Get the gene fasta by gene ID from the GTF file, and the FastaIndex
     * Before fetching the sequence, the exons are stitched together after removing overlaps
     */
    public Fasta getGeneFasta(String geneId, FastaIndex fastaIndex) {
        Feature gene = genes.get(geneId);
        if (gene == null) return null;

        // Extract, gene start and end and fetch sequences
        String sequence = getFeatureSequence(
                Stream.of(gene),
                fastaIndex
        );

        // Reverse complement if on the negative strand
        if (gene.strand.equals("-")) {
            sequence = reverseComplement(sequence);
        }

        return new Fasta(geneIDs.indexOf(geneId), geneId, sequence);
    }

    /**
     * Get the transcript fasta by transcript ID from the GTF file, and the FastaIndex
     */
    public Fasta getTranscriptFasta(String transcriptId, FastaIndex fastaIndex) {
        Feature transcript = findFeatureById(transcriptId, "transcript");
        if (transcript == null) return null;

        // Extract, stitch, and fetch sequences for all exons of the transcript
        String sequence = getFeatureSequence(
                transcript.children.stream().filter(f -> f.type.equals("exon")),
                fastaIndex
        );

        // Reverse complement if on the negative strand
        if (transcript.strand.equals("-")) {
            sequence = reverseComplement(sequence);
        }

        return new Fasta(transcriptIDs.indexOf(transcriptId), transcriptId, sequence);
    }

    /**
     * Get the exon fasta by exon ID from the GTF file, and the FastaIndex
     */
    public Fasta getExonFasta(String exonId, FastaIndex fastaIndex) {
        Feature exon = findFeatureById(exonId, "exon");
        if (exon == null) return null;

        // Fetch the sequence directly
        String sequence = exon.getSequence(fastaIndex);
        if (exon.strand.equals("-")) {
            sequence = reverseComplement(sequence);
        }

        return new Fasta(-1, exonId, sequence);
    }

    /**
     * Get the CDS fasta by CDS ID from the GTF file, and the FastaIndex
     */
    public Fasta getCDSFasta(String cdsId, FastaIndex fastaIndex) {
        Feature cds = findFeatureById(cdsId, "CDS");
        if (cds == null) return null;

        // Fetch the sequence directly
        String sequence = cds.getSequence(fastaIndex);
        if (cds.strand.equals("-")) {
            sequence = reverseComplement(sequence);
        }
        return new Fasta(-1, cdsId, sequence);
    }

    /**
     * Find a feature by ID and type
     */
    private Feature findFeatureById(String id, String type) {
        return switch (type.toLowerCase()) {
            case "gene" -> genes.get(id);
            case "transcript" -> genes.values().stream()
                    .flatMap(g -> g.children.stream())
                    .filter(t -> t.id.equals(id))
                    .findFirst().orElse(null);
            case "exon" -> genes.values().stream()
                    .flatMap(g -> g.children.stream())
                    .flatMap(t -> t.children.stream())
                    .filter(f -> f.id.equals(id) && f.type.equals("exon"))
                    .findFirst().orElse(null);
            case "cds" -> genes.values().stream()
                    .flatMap(g -> g.children.stream())
                    .flatMap(t -> t.children.stream())
                    .filter(f -> f.id.equals(id) && f.type.equals("CDS"))
                    .findFirst().orElse(null);
            default -> null;
        };
    }

    /**
     * Get the sequence of a feature by stitching together the sequences of its children
     */
    private String getFeatureSequence(Stream<Feature> features, FastaIndex fastaIndex) {
        Feature[] stitchedFeatures = stitchOverlappingFeatures(
                features.toArray(Feature[]::new)
        );

        return Arrays.stream(stitchedFeatures)
                .map(f -> f.getSequence(fastaIndex))
                .reduce("", String::concat);
    }

    /**
     * Stitch together overlapping features
     */
    private Feature[] stitchOverlappingFeatures(Feature[] features) {
        if (features.length == 0) return new Feature[0];

        List<Feature> stitched = new ArrayList<>();
        Arrays.sort(features, Comparator.comparingInt(f -> f.start));

        Feature current = features[0];
        for (int i = 1; i < features.length; i++) {
            Feature next = features[i];
            if (current.end >= next.start) {
                current.end = Math.max(current.end, next.end);
            } else {
                stitched.add(current);
                current = next;
            }
        }
        stitched.add(current);

        return stitched.toArray(new Feature[0]);
    }

    /**
     * Get all features of a specific type on a specific sequence (chromosome)
     */
    public Feature[] getFeatures(String sequenceName, String featureType) {
        // Get all features of a specific type on a specific sequence (chromosome)
        return switch (featureType.toLowerCase()) {
            case "gene" -> genes.values().stream()
                    .filter(g -> g.chrom.equals(sequenceName))
                    .toArray(Feature[]::new);
            case "transcript" -> genes.values().stream()
                    .filter(g -> g.chrom.equals(sequenceName))
                    .flatMap(g -> g.children.stream())
                    .toArray(Feature[]::new);
            case "exon" -> genes.values().stream()
                    .filter(g -> g.chrom.equals(sequenceName))
                    .flatMap(g -> g.children.stream())
                    .flatMap(t -> t.children.stream())
                    .filter(f -> f.type.equals("exon"))
                    .toArray(Feature[]::new);
            case "cds" -> genes.values().stream()
                    .filter(g -> g.chrom.equals(sequenceName))
                    .flatMap(g -> g.children.stream())
                    .flatMap(t -> t.children.stream())
                    .filter(f -> f.type.equals("CDS"))
                    .toArray(Feature[]::new);
            default -> new Feature[0];
        };
    }

    /**
     * Feature class representing a GTF feature
     */
    public static class Feature {
        String type;        // Type of feature (gene, transcript, exon, CDS, etc.)
        String chrom;       // Chromosome
        int start;          // Start position
        int end;            // End position
        String strand;      // Strand (+ or -)
        String id;          // Unique identifier (e.g., gene_id, transcript_id)
        Feature parent;     // Parent feature (e.g., Gene for a Transcript)
        List<Feature> children = new ArrayList<>(); // Children features (e.g., Transcripts for a Gene)

        public Feature(String type, String chrom, int start, int end, String strand, String id) {
            this.type = type;
            this.chrom = chrom;
            this.start = start;
            this.end = end;
            this.strand = strand;
            this.id = id;
        }

        public void addChild(Feature child) {
            children.add(child);
        }

        public void sortChildrenByStart() {
            children.sort(Comparator.comparingInt(f -> f.start));
        }

        public String getSequence(FastaIndex fastaIndex) {
            return fastaIndex.getSequence(chrom, start - 1, end - start + 1);
        }

        public String getId() {
            return id;
        }

        public String getChrom() {
            return chrom;
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

        @Override
        public String toString() {
            return String.format("%s [%s:%d-%d, %s, id=%s, parent=%s, children=%d]",
                    type, chrom, start, end, strand, id, parent != null ? parent.id : "null", children.size());
        }
    }

    /**
     * Iterater over the genes
     */
    public Iterable<Feature> genes() {
        return genes.values();
    }

    /**
     * Iterater over the transcripts
     */
    public Iterable<Feature> transcripts() {
        List<Feature> transcripts = new ArrayList<>();
        for (Feature gene : genes.values()) {
            transcripts.addAll(gene.children);
        }
        return transcripts;
    }
}
// EOF