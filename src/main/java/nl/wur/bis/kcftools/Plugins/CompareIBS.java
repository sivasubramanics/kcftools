package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.Fasta;
import nl.wur.bis.kcftools.Data.FastaIndex;
import nl.wur.bis.kcftools.Data.KMC;
import nl.wur.bis.kcftools.Data.Kmer;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import picocli.CommandLine.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

@Command(name = "compareIBS", description = "Compare IBS windows between two mapping and build a all vs all matrix")
public class CompareIBS implements Callable<Integer>, Runnable {
    @Option(names = {"--refOne"}, description = "Reference one file name", required = true)
    private String refOne;
    @Option(names = {"--refTwo"}, description = "Reference two file name", required = true)
    private String refTwo;
    @Option(names = {"--kcfOne"}, description = "findIBS summary output for reference one", required = true)
    private String ibsSummaryFileOne;
    @Option(names = {"--kcfTwo"}, description = "findIBS summary output for reference two", required = true)
    private String ibsSummaryFileTwo;
    @Option(names = {"--kmc"}, description = "KMC file prefix", required = true)
    private String kmcPrefix;
    @Option(names = {"--output"}, description = "Output file name", required = true)
    private String outFile;
    @Option(names = {"-t", "--threads"}, description = "Number of threads to use", required = false)
    private int nThreads = 2;

    public CompareIBS() {
    }

    public CompareIBS(String refOne, String refTwo, String ibsSummaryFileOne, String ibsSummaryFileTwo, String kmcPrefix, String outFile, int nThreads) {
        this.refOne = refOne;
        this.refTwo = refTwo;
        this.ibsSummaryFileOne = ibsSummaryFileOne;
        this.ibsSummaryFileTwo = ibsSummaryFileTwo;
        this.kmcPrefix = kmcPrefix;
        this.outFile = outFile;
        this.nThreads = nThreads;
    }

    @Override
    public Integer call() throws Exception {
        compareIBS();
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

    private void compareIBS() throws IOException {
        try {
            HashMap<String, ArrayList<String[]>> refOneIBS = readIBSsummaryFile(ibsSummaryFileOne);
            HashMap<String, ArrayList<String[]>> refTwoIBS = readIBSsummaryFile(ibsSummaryFileTwo);
            FastaIndex indexOne = new FastaIndex(refOne);
            FastaIndex indexTwo = new FastaIndex(refTwo);
            KMC kmc = new KMC(kmcPrefix);
            BufferedWriter writer = new BufferedWriter(new java.io.FileWriter(outFile));
            for (String sample: refOneIBS.keySet()) {
                if (!refTwoIBS.containsKey(sample)) {
                    continue;
                }
                ArrayList<String[]> refOneList = refOneIBS.get(sample);
                ArrayList<String[]> refTwoList = refTwoIBS.get(sample);
                ExecutorService executor = Executors.newFixedThreadPool(nThreads);
                ExecutorCompletionService<Void> completionService = new ExecutorCompletionService<>(executor);
                List<String[]> statsList = new ArrayList<>();
                for (String[] refOneFields: refOneList) {
                    for (String[] refTwoFields: refTwoList) {
                        completionService.submit(() -> {
                                    String[] compareStats = getCommonKmers(refOneFields, indexOne, refTwoFields, indexTwo, kmc);
                                    synchronized (statsList) {
                                        statsList.add(compareStats);
                                    }
                                    return null;
                                });
                    }
                }
                for (int i = 0; i < refOneList.size() * refTwoList.size(); i++) {
                    try {
                        completionService.take().get();
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }
                }
                for (String[] stats : statsList) {
                    String line = String.join("\t", stats);
                    writer.write(line + "\n");
                }
            }
            writer.close();
            kmc.close();
            indexOne.close();
            indexTwo.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private String[] getCommonKmers(String[] refOneFields, FastaIndex indexOne, String[] refTwoFields, FastaIndex indexTwo, KMC kmc) {
        int[] stats = new int[]{0, 0, 0, 0};
        // get fasta for refOne for the start and end
        if (!indexOne.containsSequence(refOneFields[2])) {
            HelperFunctions.log("error", "CompareIBS", "Sequence " + refOneFields[2] + " not found in reference one");
        }
        HelperFunctions.log("info", "CompareIBS", "Comparing " + refOneFields[2] + " and " + refTwoFields[2]);
        Fasta fastaOne = new Fasta(Integer.parseInt(refOneFields[0]) , refOneFields[2], indexOne.getSequence(refOneFields[2], Integer.parseInt(refOneFields[3]), Integer.parseInt(refOneFields[5])));
        List<Kmer> kmersOne = fastaOne.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);
        stats[0] = kmersOne.size();
        Fasta fastaTwo = new Fasta(Integer.parseInt(refTwoFields[0]) , refTwoFields[2], indexTwo.getSequence(refTwoFields[2], Integer.parseInt(refTwoFields[3]), Integer.parseInt(refTwoFields[5])));
        List<Kmer> kmersTwo = fastaTwo.getKmersList(kmc.getKmerLength(), kmc.getPrefixLength(), false);
        stats[1] = kmersTwo.size();
        HelperFunctions.log("info", "CompareIBS", "Total kmers " + kmersOne.size() + " and " + kmersTwo.size());
        List<Kmer> commonKmers = getUniqueCommonKmersSorted(kmersOne, kmersTwo);
        HelperFunctions.log("info", "CompareIBS", "Common kmers " + commonKmers.size());
        stats[2] = commonKmers.size();
        for (Kmer kmer: commonKmers) {
            if (kmc.isExist(kmer)) {
                stats[3]++;
            }
        }
        return new String[]{refOneFields[2], refOneFields[3], refOneFields[4], refTwoFields[2], refTwoFields[3], refTwoFields[4], Integer.toString(stats[0]), Integer.toString(stats[1]), Integer.toString(stats[2]), Integer.toString(stats[3])};
    }

    public static List<Kmer> getUniqueCommonKmersSorted(List<Kmer> list1, List<Kmer> list2) {
        // Determine the smaller list for HashSet conversion to minimize memory usage
        List<Kmer> smallerList = (list1.size() <= list2.size()) ? list1 : list2;
        List<Kmer> largerList = (list1.size() > list2.size()) ? list1 : list2;

        // Use a HashSet for fast lookup and avoid duplicates
        HashSet<Kmer> set = new HashSet<>(smallerList);
        Set<Kmer> commonElements = new HashSet<>();

        // Iterate over the larger list and check for common elements
        for (Kmer kmer : largerList) {
            if (set.contains(kmer)) {
                commonElements.add(kmer); // Add to the result set
            }
        }

        // Convert the set of common elements to a list
        return new ArrayList<>(commonElements);
    }

    private HashMap<String, ArrayList<String[]>> readIBSsummaryFile(String ibsSummaryFile) {
        HashMap<String, ArrayList<String[]>> summaryIBS = new HashMap<String, ArrayList<String[]>>();
        try (BufferedReader reader = new BufferedReader(new FileReader(ibsSummaryFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields[0].equals("Block")) {
                    continue;
                }
                String sample = fields[1];
                if (summaryIBS.containsKey(sample)) {
                    summaryIBS.get(sample).add(fields);
                } else {
                    // create new list for the sample as value and add fields as first element
                    ArrayList<String[]> list = new ArrayList<String[]>();
                    list.add(fields);
                    summaryIBS.put(sample, list);
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return summaryIBS;
    }

}
