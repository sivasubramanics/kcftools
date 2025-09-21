package nl.wur.bis.kcftools.Plugins;

import nl.wur.bis.kcftools.Data.KCFHeader;
import nl.wur.bis.kcftools.Data.KCFReader;
import nl.wur.bis.kcftools.Data.KCFWriter;
import nl.wur.bis.kcftools.Data.Window;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine.*;

import java.io.IOException;
import java.util.concurrent.Callable;

@Command(name = "scoreRecalc", description = "Recalculate scores in a KCF file")
public class ScoreRecalc implements Callable<Integer>, Runnable {
    @Option(names = {"-i", "--input"}, description = "Input KCF file", required = true)
    private String inFile;

    @Option(names = {"-o", "--output"}, description = "Output KCF file", required = true)
    private String outFile;

    // inner kmer distance weight
    @Option(names = {"--wi"}, description = "Inner kmer distance weight [0.3]", required = false)
    private double innerDistanceWeight = 0.3;
    // tail kmer distance weight
    @Option(names = {"--wt"}, description = "Tail kmer distance weight [0.3]", required = false)
    private double tailDistanceWeight = 0.3;
    // kmer ratio weight
    @Option(names = {"--wr"}, description = "Kmer ratio weight [0.4]", required = false)
    private double kmerRatioWeight = 0.4;

    private static final String CLASS_NAME = ScoreRecalc.class.getSimpleName();

    @Override
    public Integer call() throws Exception {
        double[] weights = getWeights();
        recalculateScores(weights);
        return 0;
    }

    @Override
    public void run() {
        try {
            call();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void recalculateScores(double[] weights) {
        try (KCFReader reader = new KCFReader(inFile);
             KCFWriter writer = new KCFWriter(outFile)) {
            KCFHeader header = reader.getHeader();
            header.setWeightInnerDist(innerDistanceWeight);
            header.setWeightTailDist(tailDistanceWeight);
            header.setWeightKmerRatio(kmerRatioWeight);
            writer.writeHeader(header);
            for (Window window : reader) {
                window.recalcScore(weights);
                writer.writeWindow(window);
            }
            Logger.info(CLASS_NAME, "Recalculated scores and wrote to " + outFile);
        } catch (IOException e) {
            throw new RuntimeException("Error processing KCF files", e);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private double[] getWeights() {
        return new double[]{innerDistanceWeight, tailDistanceWeight, kmerRatioWeight};
    }
}
