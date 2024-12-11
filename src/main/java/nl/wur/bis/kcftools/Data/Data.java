package nl.wur.bis.kcftools.Data;

public class Data{
    int kmerDistance;
    int observedKmers;
    int variations;
    float score;
    int ibs;

    public Data(int observedKmers, int variations, int kmerDistance, int totalKmers){
        this.observedKmers = observedKmers;
        this.variations = variations;
        this.kmerDistance = kmerDistance;
        this.score = computeScore(totalKmers);
        this.ibs = -1;
    }

    public Data() {
        this.ibs = -1;
    }

    public float computeScore(int totalKmers){
        if ((observedKmers == 0) || (totalKmers == 0)){
            return 0;
        }
//        return (float) (observedKmers - variations) * 100 / totalKmers;
//        var$score <- (var$observed_kmers / var$total_kmers) * (1 - (var$kmer_distance / ((var$end-var$start+1)))) * 100
        return (float) observedKmers / totalKmers * (1 - (float) kmerDistance / (float) (totalKmers)) * 100;

    }

    public float getScore(){
        return score;
    }

    @Override
    public String toString(){
        if (ibs == -1){
            return "N" + ":" + variations + ":" + observedKmers + ":" + kmerDistance + ":" + String.format("%.2f", score);
        }
        return ibs + ":" + variations + ":" + observedKmers + ":" + kmerDistance + ":" + String.format("%.2f", score);
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
        return observedKmers + "\t" + variations + "\t" + kmerDistance;
    }
}
// EOF