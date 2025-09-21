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
    int tailDistance;
    int innerDistance;
    int observedKmers;
    int variations;
    double score;
    int ibs;
    int leftDistance;
    int rightDistance;

    private final String CLASS_NAME = this.getClass().getSimpleName();

    public Data(int observedKmers,
                int variations,
                int innerDistance,
                int leftDistance,
                int rightDistance,
                int totalKmers,
                int effLength,
                double[] weights) {
        this.observedKmers = observedKmers;
        this.variations = variations;
        this.leftDistance = leftDistance;
        this.rightDistance = rightDistance;
        this.innerDistance = innerDistance;
        this.score = computeScore(totalKmers, effLength, weights);
        this.ibs = -1;
    }

    public Data(int observedKmers,
                int variations,
                int innerDistance,
                int leftDistance,
                int rightDistance,
                int totalKmers,
                int effLength,
                String ibs, double[] weights) {
        this.observedKmers = observedKmers;
        this.variations = variations;
        this.leftDistance = leftDistance;
        this.rightDistance = rightDistance;
        this.innerDistance = innerDistance;
        this.score = computeScore(totalKmers, effLength, weights);
        if (ibs.equals("N")){
            this.ibs = -1;
        } else {
            this.ibs = Integer.parseInt(ibs);
        }
    }

    /***
     * Compute the score of the window
     */
    public double computeScore(int totalKmers, int effLength, double[] weights){
        if ((observedKmers == 0) || (totalKmers == 0) || (effLength == 0)){
            return 0;
        }
//        old method of scoring
//        return (float) (observedKmers - variations) * 100 / totalKmers;
        return ((weights[2] * ((double) observedKmers / totalKmers))
                + (weights[0] * (1.0f - ((double) innerDistance / effLength)))
                + (weights[1] * (1.0f - ((double) tailDistance / effLength)))) * 100.0f;
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
    public String toString(){
        // ibs:variations:observedKmers:innerDistance:leftDistane:rightDistance:score
        if (ibs == -1){
            return "N" + ":" + variations + ":" + observedKmers + ":" + innerDistance + ":" + leftDistance + ":" + rightDistance + ":" + String.format("%.2f", score);
        }
        return ibs + ":" + variations + ":" + observedKmers + ":" + innerDistance + ":" + leftDistance + ":" + rightDistance + ":" + String.format("%.2f", score);
    }

    /***
     * Get TSV formated string representation of the data (to write in native IBSpy table format)
     */
    public String toTSV() {
        int distance = getInnerDistance() + getTailDistance();
        double score = getScore();
        return getObservedKmers() + "\t" + getVariations() + "\t" + distance + "\t" + String.format("%.2f", score);
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
}
// EOF