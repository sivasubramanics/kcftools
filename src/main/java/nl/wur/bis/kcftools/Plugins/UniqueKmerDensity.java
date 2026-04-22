package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.*;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;


@Command(name = "uniqueKmerDensity", description = "Calculate the density of unique kmers in a fasta sequences")
public class UniqueKmerDensity implements Callable<Integer>, Runnable {
    @Option(names = {"-f", "--fasta"}, description = "Input fasta file", required = true)
    private String fastaFile;

    @Option(names = {"-k", "--kmc"}, description = "KMC database prefix (created from the same fasta file)", required = true)
    private String kmcPrefix;

    @Option(names = {"-w", "--window"}, description = "Window size (default: 100000)", defaultValue = "100000")
    private int windowSize;

    @Option(names = {"-m", "--memory"}, description = "Load KMC database into memory (default: false)", defaultValue = "false")
    private boolean loadKmcInMemory;

    @Option(names = {"-o", "--output"}, description = "Output file name (default: unique_kmer_density.txt)", defaultValue = "unique_kmer_density.txt")
    private String outFile;

    @Option(names = {"-t", "--threads"}, description = "Number of threads to use (default: auto-detected)", defaultValue = "2")
    private int nThreads;

    private static final String CLASS_NAME = UniqueKmerDensity.class.getSimpleName();
    private FastaIndex index;

    public UniqueKmerDensity() {
    }

    @Override
    public Integer call() throws Exception {
        HelperFunctions.printCommandLine(new CommandLine(this), CLASS_NAME);
        calculateUniqueKmerDensity();
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private void calculateUniqueKmerDensity() throws IOException {

        KMC kmc = new KMC(kmcPrefix, loadKmcInMemory);
        int kmerSize = kmc.getKmerLength();
        FastaIndex index = new FastaIndex(fastaFile);

        ConcurrentHashMap<String, Queue<Bin>> binMap = new ConcurrentHashMap<>();
        for (String seqName : index.getSequenceNames()) {
            int seqLength = index.getSequenceLength(seqName);
            Queue<Bin> bins = createBins(seqName, seqLength, kmerSize, windowSize);
            binMap.put(seqName, bins);
        }

        int totalBins = binMap.values().stream().mapToInt(Queue::size).sum();
        Logger.info(CLASS_NAME, "Created bins for all sequences. Total bins: " + totalBins);

        LinkedHashMap<String, List<Bin>> processedBins = new LinkedHashMap<>();
        AtomicInteger completedBins = new AtomicInteger(0);

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);

        for (Map.Entry<String, Queue<Bin>> entry : binMap.entrySet()) {
            String name = entry.getKey();
            Queue<Bin> bins = entry.getValue();
            List<Bin> processed = Collections.synchronizedList(new ArrayList<>(bins.size()));
            processedBins.put(name, processed);

            for (Bin bin : bins) {
                completionService.submit(() -> {
                    Fasta fasta = new Fasta(1, "x", index.getSequence(bin.getName(), bin.getStart(), bin.getLength()), bin.getName() + ":" + bin.getStart() + "-" + bin.getEnd());
                    Bin processedBin = processBin(bin, fasta, kmc);
                    processed.add(processedBin);
                    int completed = completedBins.incrementAndGet();
                    float progress = (float) (completed * 100) / totalBins;
                    synchronized (System.out) {
                        System.out.printf("\rProgress: %.2f%%", progress);
                    }
                    return null;
                });
            }
        }

        for (int i = 0; i < totalBins; i++) {
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

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outFile))) {
            writer.write("name\tbin_start\tbin_end\ttotal_kmers\tuniq_bin\tuniq_overall\n");
            // sort the windows based on its start position
            for (String name : processedBins.keySet()) {
                processedBins.get(name).sort(Comparator.comparingInt(Bin::getStart));
                for (Bin bin : processedBins.get(name)) {
                    writer.write(bin.toString() + "\n");
                }
            }
        } catch (Exception e) {
            Logger.error(CLASS_NAME, "Can not create output file.");
        }
    }

    private Bin processBin(Bin bin, Fasta fasta, KMC kmc) {

        int totalKmers = 0;
        int uniqueBinKmers = 0;
        int uniqueOverAllKmers = 0;

        if (fasta == null) {
            Logger.error(CLASS_NAME, "Fasta object is null for bin. " + bin.getName() + ":" + bin.getStart() + "-" + bin.getEnd());
            return bin;
        }
        List<Kmer> kmers = fasta.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);
        totalKmers = kmers.size();
        Set<Kmer> kmerSet = new HashSet<>(kmers);

        if (!kmerSet.isEmpty()) {
            for (Kmer k : kmerSet) {
                uniqueBinKmers++;
                Kmer km = new Kmer(k, kmc.isBothStrands());
                int kmerCount = kmc.getCount(km);
                if (kmerCount == 1) {
                    uniqueOverAllKmers++;
                }
            }
        }

        synchronized (bin) {
            bin.setTotalKmers(totalKmers);
            bin.setUniqueBinKmers(uniqueBinKmers);
            bin.setUniqueOverAllKmers(uniqueOverAllKmers);
        }

        return bin;
    }

    private Queue<Bin> createBins(String seqName, int seqLength, int kmerSize, int windowSize) {
        Queue<Bin> bins = new ConcurrentLinkedDeque<>();
        int lastEnd = 0;
        while (lastEnd < seqLength) {
            int start = Math.max(0, lastEnd - kmerSize + 1);
            int end = Math.min(start + windowSize, seqLength);

            if (end - start >= kmerSize) {
                bins.add(new Bin(seqName, start, end));
            }

            lastEnd = end;
        }
        return bins;
    }


    class Bin {
        private String name;
        private int start;
        private int end;
        private long totalKmers;
        private long uniqueBinKmers;
        private long uniqueOverAllKmers;

        public Bin(String name, int start, int end) {
            this.name = name;
            this.start = start;
            this.end = end;
            this.totalKmers = 0;
            this.uniqueBinKmers = 0;
            this.uniqueOverAllKmers = 0;
        }

        public String getName() {
            return name;
        }
        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public int getLength() {
            return end - start;
        }

        public long getTotalKmers() {
            return totalKmers;
        }

        public long getUniqueBinKmers() {
            return uniqueBinKmers;
        }

        public long getUniqueOverAllKmers() {
            return uniqueOverAllKmers;
        }

        public void setUniqueBinKmers(long uniqueBinKmers) {
            this.uniqueBinKmers = uniqueBinKmers;
        }

        public void setUniqueOverAllKmers(long uniqueOverAllKmers) {
            this.uniqueOverAllKmers = uniqueOverAllKmers;
        }

        public void setTotalKmers (long totalKmers) {
            this.totalKmers = totalKmers;
        }

        @Override
        public String toString() {
            return String.format("%s\t%d\t%d\t%d\t%d\t%d",
                    name, start, end, totalKmers, uniqueBinKmers, uniqueOverAllKmers);
        }
    }

}


