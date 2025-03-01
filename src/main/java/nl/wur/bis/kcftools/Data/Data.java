package nl.wur.bis.kcftools.Data;

public class Data{
    int tailDistance;
    int innerDistance;
    int observedKmers;
    int variations;
    double score;
    int ibs;

    public Data(int observedKmers, int variations, int innerDistance, int tailDistance, int totalKmers, int effLength, double[] weights) {
        this.observedKmers = observedKmers;
        this.variations = variations;
        this.innerDistance = innerDistance;
        this.tailDistance = tailDistance;
        this.score = computeScore(totalKmers, effLength, weights);
        this.ibs = -1;
    }

    public Data(int observedKmers, int variations, int innerDistance, int tailDistance, int totalKmers, int effLength, String ibs, double[] weights) {
        this.observedKmers = observedKmers;
        this.variations = variations;
        this.innerDistance = innerDistance;
        this.tailDistance = tailDistance;
        this.score = computeScore(totalKmers, effLength, weights);
        if (ibs.equals("N")){
            this.ibs = -1;
        } else {
            this.ibs = Integer.parseInt(ibs);
        }
    }

    public double computeScore(int totalKmers, int effLength, double[] weights){
        if ((observedKmers == 0) || (totalKmers == 0) || (effLength == 0)){
            return 0;
        }
//        return (float) (observedKmers - variations) * 100 / totalKmers;
//        var$score <- (var$observed_kmers / var$total_kmers) * (1 - (var$kmer_distance / ((var$end-var$start+1)))) * 100
//        return (float) observedKmers / totalKmers * (1 - (float) kmerDistance / (float) (totalKmers)) * 100;
        return ((weights[2] * ((double) observedKmers / totalKmers))
                + (weights[0] * (1.0f - ((double) innerDistance / effLength)))
                + (weights[1] * (1.0f - ((double) tailDistance / effLength)))) * 100.0f;
    }

    public double getScore(){
        return score;
    }

    @Override
    public String toString(){
        // ibs:variations:observedKmers:innerDistance:tailDistance:score
        if (ibs == -1){
            return "N" + ":" + variations + ":" + observedKmers + ":" + innerDistance + ":" + tailDistance + ":" + String.format("%.2f", score);
        }
        return ibs + ":" + variations + ":" + observedKmers + ":" + innerDistance + ":" + tailDistance + ":" + String.format("%.2f", score);
    }

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

    public String toTSV() {
        return observedKmers + "\t" + variations + "\t" + innerDistance;
    }

    public int getInnerDistance() {
        return innerDistance;
    }

    public int getTailDistance() {
        return tailDistance;
    }
}
// EOF