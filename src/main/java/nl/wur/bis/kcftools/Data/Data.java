package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Logger;
/**
 * Data class to store the data for each window
 * observedKmers: number of kmers observed in the window
 * variations: number of variations in the window (instances of stretch of bases that are not part of observedKmers)
 * innerDistance: cumulative number of bases that are not part of the observedKmers, which are between the left and right observedKmers
 * leftDistance: number of bases that are not part of the observedKmers, which are to the left side of window
 * rightDistance: number of bases that are not part of the observedKmers, which are to the right side of window
 * score: identity score of the window
 * ibs: N/1 presense or absence of IBS
 * tailDistance: number of bases that are not part of the observedKmers, which are to the left and right side of window
 */
public class Data{
    int innerDistance;
    int observedKmers;
    int variations;
    double score;
    int ibs;
    int leftDistance;
    int rightDistance;
    double meanKmerCount;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public Data(int observedKmers,
                int variations,
                int innerDistance,
                int leftDistance,
                int rightDistance,
                long kmerCount,
                int totalKmers,
                int effLength,
                double[] weights) {
        this(observedKmers, variations, innerDistance, leftDistance, rightDistance,
                kmerCount, totalKmers, effLength, "-1", weights);
    }

    public Data(int observedKmers,
                int variations,
                int innerDistance,
                int leftDistance,
                int rightDistance,
                long kmerCount,
                int totalKmers,
                int effLength,
                String ibs,
                double[] weights) {

        this.observedKmers = observedKmers;
        this.variations = variations;
        this.leftDistance = leftDistance;
        this.rightDistance = rightDistance;
        this.innerDistance = innerDistance;

        this.meanKmerCount = (kmerCount > 0)
                ? (double) kmerCount / observedKmers
                : 0.00;

        this.score = computeScore(totalKmers, effLength, weights);

        if ("N".equals(ibs)) {
            this.ibs = -1;
        } else {
            this.ibs = Integer.parseInt(ibs);
        }
    }

    public void update(int observedKmers,
                       int variations,
                       int innerDistance,
                       int leftDistance,
                       int rightDistance,
                       long kmerCount,
                       String ibs,
                       int totalKmers,
                       int effLength,
                       double[] weights) {

        this.observedKmers = observedKmers;
        this.variations = variations;
        this.innerDistance = innerDistance;
        this.leftDistance = leftDistance;
        this.rightDistance = rightDistance;

        this.meanKmerCount = (observedKmers > 0) ? (double) kmerCount / observedKmers : 0.00;
        this.score = computeScore(totalKmers, effLength, weights);
        this.ibs = "N".equals(ibs) ? -1 : Integer.parseInt(ibs);
    }

    /***
     * Compute the score of the window
     */
    public double computeScore(int totalKmers, int effLength, double[] weights){
        if ((getObservedKmers() == 0) || (totalKmers == 0) || (effLength == 0)){
            return 0;
        }
//        old method of scoring
//        return (float) (observedKmers - variations) * 100 / totalKmers;
        if (weights[0] + weights[1] + weights[2] != 1.0){
            Logger.error(CLASS_NAME, "Weights should sum to 1.0");
        }
        return ((weights[2] * ((double) getObservedKmers() / totalKmers))
                + (weights[0] * (1.0f - ((double) getInnerDistance() / effLength)))
                + (weights[1] * (1.0f - ((double) getTailDistance() / effLength)))) * 100.0f;
    }

    /***
     * Get the score of the window
     */
    public double getScore(){
        return score;
    }

    /***
     * Get string representation of the data (for KCF writing purpose)
     */
    @Override
    public String toString() {
        String ibsValue = (ibs == -1) ? "N" : String.valueOf(ibs);
        return String.join(":",
                ibsValue,
                String.valueOf(getVariations()),
                String.valueOf(getObservedKmers()),
                String.valueOf(getInnerDistance()),
                String.valueOf(getLeftDistance()),
                String.valueOf(getRightDistance()),
                String.format("%.2f", getMeanKmerCount()),
                String.format("%.2f", getScore())
        );
    }

    /***
     * Get TSV formated string representation of the data (to write in native IBSpy table format)
     */
    public String toTSV() {
        return String.join("\t",
                String.valueOf(getObservedKmers()),
                String.valueOf(getVariations()),
                String.valueOf(getInnerDistance() + getTailDistance()),
                String.format("%.2f", getMeanKmerCount()),
                String.format("%.2f", getScore())
        );
    }

    /***
     * Setter and getter functions follows
     */
    public void setIBS(int ibs) {
        this.ibs = ibs;
    }

    public int getIBS() {
        return ibs;
    }

    public int getVariations() {
        return variations;
    }

    public int getObservedKmers() {
        return observedKmers;
    }

    public int getInnerDistance() {
        return innerDistance;
    }

    public int getTailDistance() {
        return getLeftDistance() + getRightDistance();
    }

    public int getLeftDistance() {
        return leftDistance;
    }

    public int getRightDistance() {
        return rightDistance;
    }

    public double getMeanKmerCount() {
        return meanKmerCount;
    }
}
// EOF