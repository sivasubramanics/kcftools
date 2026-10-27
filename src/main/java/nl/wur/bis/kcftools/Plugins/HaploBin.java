package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

@Command(name = "haplobin", description = "Find the Haplotype bins from hetrotic genomes", sortOptions = false)
public class HaploBin implements Runnable, Callable<Integer> {
    @Option(names = {"-a", "--kmc_a"}, description = "kmc DB prefix of genome A.", required = true)
    private String kmcDBPrefixA;

    @Option(names = {"-b", "--kmc_b"}, description = "kmc DB prefix of genome B.", required = true)
    private String kmcDBPrefixB;

    @Option(names = {"-x", "--kmc_x"}, description = "kmc DB prefix of genome X (F1).", required = true)
    private String kmcDBPrefixX;

    @Option(names = {"-X", "--genome_x"}, description = "Path to the genome X (F1) fasta file.", required = true)
    private String genomeX;

    @Option(names = {"-o", "--output"}, description = "Path to the output file.", required = true)
    private String outputFile;

    @Option(names = {"-w", "--window_size"}, description = "Window size for k-mer analysis.", required = false, defaultValue = "50000")
    private int windowSize;

    @Option(names = {"-s", "--step_size"}, description = "Step size for sliding window analysis. If set to 0, non-overlapping windows will be used.", required = false, defaultValue = "0")
    private int stepSize;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use for processing.", required = false, defaultValue = "2")
    private int nThreads;

    private final String CLASS_NAME = HaploBin.class.getSimpleName();

    int kmerSize;
    int prefixSize;
    FastaIndex index;

    public HaploBin() {

    }

    @Override
    public void run() {
        // Implementation for running the HaploBin tool
        try {
            call();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public Integer call() throws Exception {
        // Implementation for calling the HaploBin tool
        HelperFunctions.printCommandLine(new CommandLine(this), this.getClass().getSimpleName());
        getHaploBins();
        return 0;
    }

    private  void getHaploBins() throws Exception {
        try {
            KMC kmcA = new KMC(kmcDBPrefixA);
            KMC kmcB = new KMC(kmcDBPrefixB);
            KMC kmcX = new KMC(kmcDBPrefixX);
            if (kmcA.getKmerLength() != kmcB.getKmerLength() || kmcA.getKmerLength() != kmcX.getKmerLength()) {
                Logger.error(CLASS_NAME, "K-mer lengths of the provided KMC databases do not match. Please ensure that all KMC databases were generated with the same k-mer size.");
            }
            if (kmcA.getPrefixLength() != kmcB.getPrefixLength() || kmcA.getPrefixLength() != kmcX.getPrefixLength()) {
                Logger.error(CLASS_NAME, "KMC prefix lengths of the provided KMC databases do not match. Please ensure that all KMC databases were generated with the same prefix length.");
            }
            kmerSize = kmcA.getKmerLength();
            prefixSize = kmcA.getPrefixLength();
            index = new FastaIndex(genomeX);
            ConcurrentHashMap<String, Queue<Region>> regionsMap = new ConcurrentHashMap<>();
            for (String name : index.getSequenceNames()) {
                Queue<Region> regions = getRegions(name);
                regionsMap.put(name, regions);
            }
            int totalRegions = regionsMap.values().stream().mapToInt(Queue::size).sum();
            Logger.info(CLASS_NAME, "Total regions to process: " + totalRegions);

            LinkedHashMap<String, List<Region>> processedRegions = new LinkedHashMap<>();
            AtomicInteger completedRegions = new AtomicInteger(0);

            ExecutorService executor = Executors.newFixedThreadPool(nThreads);
            ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);

            for (Map.Entry<String, Queue<Region>> entry : regionsMap.entrySet()) {
                String name = entry.getKey();
                Queue<Region> regions = entry.getValue();
                List<Region> processed = Collections.synchronizedList(new ArrayList<>(regions.size()));
                processedRegions.put(name, processed);

                for (Region region : regions) {
                    completionService.submit(() -> {
                        Region processedRegion = processRegion(region, kmcA, kmcB, kmcX, index);
                        processed.add(processedRegion);
                        int completed = completedRegions.incrementAndGet();
                        float progress = (float) (completed * 100) / totalRegions;
                        synchronized (System.out) {
                            System.out.printf("\rProgress: %.2f%%", progress);
                        }
                        return null;
                    });
                }
            }
            for (int i = 0; i < totalRegions; i++) {
                try {
                    completionService.take().get();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            System.out.print("\r");
            for (int i = 0; i < 100; i++) {
                System.out.print(" ");
            }
            System.out.print("\r");
            System.out.flush();

            // sort the regions based on its start position
            for (String name : processedRegions.keySet()) {
                processedRegions.get(name).sort(Comparator.comparingInt(Region::getStart));
            }

            // write the processed regions to the output file

            try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputFile))) {
                writer.write("RegionName\tSequenceName\tStart\tEnd\tTotalUniqueKmers\tKmersX\tKmersA\tKmersB\n");
                for (List<Region> regions : processedRegions.values()) {
                    for (Region region : regions) {
                        writer.write(region.toString());
                        writer.newLine();
                    }
                }
            } catch (IOException e) {
                Logger.error(CLASS_NAME, "Error writing output file: " + e.getMessage());
            }
            index.close();
            kmcA.close();
            kmcB.close();
            kmcX.close();
            HelperFunctions.printMaxMemoryUsage();
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new Exception("Error initializing KMC databases: " + e.getMessage());
        }

    }

    private Region processRegion(Region region, KMC kmcA, KMC kmcB, KMC kmcX, FastaIndex index) {
        // Implementation for processing a single region
        // This method will analyze the k-mer counts in the specified region for genomes A, B, and X
        // and update the region's k-mer counts accordingly.
        Fasta fasta = new Fasta(1, region.getName(), index.getSequence(region.getSequenceName()).substring(region.getStart(), region.getEnd()));
        int totalUniqKmers = 0;
        int kmersA = 0;
        int kmersB = 0;
        int kmersX = 0;
        List<Kmer> kmers = fasta.getKmersList(kmerSize, prefixSize, false);
        if (!kmers.isEmpty()) {
            for (Kmer k : kmers) {
                // 1. Must be unique in X
                int countX = kmcX.getCount(k);
                if (countX != 1) {
                    continue;
                }
                region.addTotalUniqKmers(1);

                int countA = kmcA.getCount(k);
                int countB = kmcB.getCount(k);

                // 2. Must be absent in one parent and present in the other
                // (You can also adjust thresholds here, e.g., countA > 1 to avoid parental errors)
                boolean inA = (countA > 0);
                boolean inB = (countB > 0);

                if (inA && inB) {
                    continue; // Shared between both parents, not informative
                }
                if (!inA && !inB) {
                    continue; // Present in X but neither parent (sequencing error/artifact) -> SKIP
                }

                // Now we know it is strictly informative (either A-exclusive or B-exclusive)
                region.addKmersX(1);

                if (inA) {
                    region.addKmersA(1);
                } else if (inB) {
                    region.addKmersB(1);
                }
            }
        }
        return region;
    }


    private Queue<Region> getRegions(String sequenceName) {
        Queue<Region> regions = new ConcurrentLinkedDeque<>();
        int sequenceLength = index.getSequenceLength(sequenceName);
        int lastPos = 0;

        if (stepSize > 0) {
            // Sliding window mode
            while (lastPos < sequenceLength) {
                int start = lastPos;
                int end = Math.min(start + windowSize, sequenceLength);

                if (end - start >= kmerSize) {
                    regions.add(new Region(sequenceName + "_" + start, sequenceName, start, end));
                }

                lastPos += stepSize; // slide by stepSize
            }
        } else {
            // Old tiling mode (non-overlapping but with kmer overlap at boundary)
            int lastEnd = 0;
            while (lastEnd < sequenceLength) {
                int start = Math.max(0, lastEnd - kmerSize + 1);
                int end = Math.min(start + windowSize, sequenceLength);

                if (end - start >= kmerSize) {
                    regions.add(new Region(sequenceName + "_" + start, sequenceName, start, end));
                }

                lastEnd = end;
            }
        }
        return regions;
    }

    static class Region {
        String name;
        String sequenceName;
        int start;
        int end;
        int totalUniqKmers;
        int kmersX;
        int kmersB;
        int kmersA;


        public Region(String name, String sequenceName, int start, int end) {
            this.name = name;
            this.sequenceName = sequenceName;
            this.start = start;
            this.end = end;
            this.totalUniqKmers = 0;
            this.kmersX = 0;
            this.kmersB = 0;
            this.kmersA = 0;
        }

        public int getKmersA() {
            return kmersA;
        }

        public int getKmersB() {
            return kmersB;
        }

        public int getKmersX() {
            return kmersX;
        }
        public int getTotalUniqKmers() {
            return totalUniqKmers;
        }
        public void addKmersA(int count) {
            this.kmersA += count;
        }
        public void addKmersB(int count) {
            this.kmersB += count;
        }
        public void addKmersX(int count) {
            this.kmersX += count;
        }
        public void addTotalUniqKmers(int count) {
            this.totalUniqKmers += count;
        }
        @Override
        public String toString() {
            return name + "\t" + sequenceName + "\t" + start + "\t" + end + "\t" + totalUniqKmers + "\t" + kmersX + "\t" + kmersB + "\t" + kmersA;
        }

        public String getName() {
            return name;
        }
        public String getSequenceName() {
            return sequenceName;
        }
        public int getStart() {
            return start;
        }
        public int getEnd() {
            return end;
        }
    }

}
