package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Objects;

/***
 * Class to represent a window in the KCF file
 */
public class Window implements Comparable<Window> {
    private final String CLASSNAME = this.getClass().getSimpleName();
    String windowId;
    String sequenceName;
    int start;
    int end;
    int effLength;
    int totalKmers;
    HashMap<String, Data> data;
    int minObservedKmers;
    int maxObservedKmers;
    float meanObservedKmers;
    int minVariations;
    int maxVariations;
    float meanVariations;
    double minScore;
    double maxScore;
    double meanScore;

    public Window(String windowId, String sequenceName, int start, int end){
        this.windowId = windowId;
        this.sequenceName = sequenceName;
        this.start = start;
        this.end = end;
        this.effLength = 0;
        this.totalKmers = 0;
        this.data = new LinkedHashMap<>();
    }

    public Window(String[] fields, String[] samples, double[] weights) {
        this.sequenceName = fields[0];
        this.start = Integer.parseInt(fields[1]);
        this.end = Integer.parseInt(fields[2]);
        this.windowId = fields[3];
        this.totalKmers = Integer.parseInt(fields[4]);
        this.effLength = Integer.parseInt(getInfoFieldMap(fields[5]).get("EFFLEN"));
        this.data = new LinkedHashMap<>();
        for (int i = 7; i < fields.length; i++) {
            String sampleName = samples[i - 7];
            Data d = parseSampleData(fields[i], totalKmers, effLength, weights);
            data.put(sampleName, d);
        }
    }

    private Data parseSampleData(String field,
                                 int totalKmers,
                                 int effLength,
                                 double[] weights) {
        // sampleData: ibs:variations:observedKmers:innerDistance:tailDistance:rightDistance:kmerCount:score
        String[] sampleData = field.split(":");
        String ibs = sampleData[0];
        int variations = Integer.parseInt(sampleData[1]);
        int observedKmers = Integer.parseInt(sampleData[2]);
        int innerDistance = Integer.parseInt(sampleData[3]);
        int leftDistance = Integer.parseInt(sampleData[4]);
        int rightDistance = Integer.parseInt(sampleData[5]);
        // Ensure kmerCount is a long, not double
        long kmerCount = Math.round(Double.parseDouble(sampleData[6]) * observedKmers);
        return new Data(
                observedKmers,
                variations,
                innerDistance,
                leftDistance,
                rightDistance,
                kmerCount,
                totalKmers,
                effLength,
                ibs,
                weights
        );
    }

    /***
     * Add data to the window
     */
    public void addData(HashMap<String, Data> data){
        for (String sample : data.keySet()){
            if (this.data.containsKey(sample)){
                Logger.error(CLASSNAME, "Sample " + sample + " already exists in window " + windowId);
            } else {
                this.data.put(sample, data.get(sample));
            }
        }
    }

    /***
     * Add data to the window
     */
    public synchronized void addData(String sample,
                                     int observedKmers,
                                     int variations,
                                     int innerDistance,
                                     int leftDistance,
                                     int rightDistance,
                                     long kmerCount,
                                     String ibs,
                                     double[] weights) {

        Data d = data.computeIfAbsent(sample,
                k -> new Data(0, 0, 0, 0, 0, 0, totalKmers, effLength, weights));

        d.update(observedKmers, variations, innerDistance, leftDistance,
                rightDistance, kmerCount, ibs, totalKmers, effLength, weights);
    }

    public void recalcScore(double[] weights){
        for (Data d : data.values()){
            d.score = d.computeScore(totalKmers, effLength, weights);
        }
    }

    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        sb.append(sequenceName)
                .append("\t").append(start)
                .append("\t").append(end)
                .append("\t").append(windowId)
                .append("\t").append(totalKmers)
                .append("\t").append(getInfoField())
                .append("\t").append(getFormatField());
        for (Data d : data.values()){
            sb.append("\t").append(d.toString());
        }
        return sb.toString();
    }

    private String getInfoField() {
        calculateStats();
        return "EFFLEN=" + effLength + ";" +
                "IS=" + String.format("%.2f", minScore) + ";" +
                "XS=" + String.format("%.2f", maxScore) + ";" +
                "MS=" + String.format("%.2f", meanScore) + ";" +
                "IO=" + minObservedKmers + ";" +
                "XO=" + maxObservedKmers + ";" +
                "MO=" + String.format("%.2f", meanObservedKmers) + ";" +
                "IV=" + minVariations + ";" +
                "XV=" + maxVariations + ";" +
                "MV=" + meanVariations;
    }

    /***
     * Parse the info field to a map
     */
    private HashMap<String, String> getInfoFieldMap(String infoField) {
        HashMap<String, String> info = new LinkedHashMap<>();
        String[] fields = infoField.split(";");
        for (int i = 0; i < fields.length; i++){
            fields[i] = fields[i].replace(";", "");
        }
        for (String field : fields){
            String[] kv = field.split("=");
            info.put(kv[0], kv[1]);
        }
        return info;
    }

    private String getFormatField() {
        return "GT:VA:OB:ID:LD:RD:KD:SC";
    }

    /***
     * Calculate the stats for the window to be used in the info field
     */
    public void calculateStats(){
        minObservedKmers = Integer.MAX_VALUE;
        maxObservedKmers = Integer.MIN_VALUE;
        meanObservedKmers = 0;
        minVariations = Integer.MAX_VALUE;
        maxVariations = Integer.MIN_VALUE;
        meanVariations = 0;
        minScore = Float.MAX_VALUE;
        maxScore = Float.MIN_VALUE;
        meanScore = 0;
        for (Data d : data.values()){
            if (d.observedKmers < minObservedKmers){
                minObservedKmers = d.observedKmers;
            }
            if (d.observedKmers > maxObservedKmers){
                maxObservedKmers = d.observedKmers;
            }
            meanObservedKmers += d.observedKmers;
            if (d.variations < minVariations){
                minVariations = d.variations;
            }
            if (d.variations > maxVariations){
                maxVariations = d.variations;
            }
            meanVariations += d.variations;
            if (d.score < minScore){
                minScore = d.score;
            }
            if (d.score > maxScore){
                maxScore = d.score;
            }
            meanScore += d.score;
        }
        // round to 2 decimal places
        meanObservedKmers /= data.size();
        meanVariations /= data.size();
        meanScore /= data.size();
    }

    public int length() {
        return end - start;
    }

    public synchronized void addTotalKmers(int totalKmersCount){
        totalKmers += totalKmersCount;
    }

    public Fasta getFasta(FastaIndex index) {
        return new Fasta(1, String.valueOf(windowId), index.getSequence(sequenceName, start, length()), sequenceName + ":" + start + "-" + end);
    }

    public String getWindowId() {
        return windowId;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Window window = (Window) o;
        return start == window.start &&
                end == window.end &&
                Objects.equals(sequenceName, window.sequenceName);
    }

    @Override
    public int hashCode() {
        return Objects.hash(sequenceName, start, end);
    }

    @Override
    public int compareTo(@NotNull Window o) {
        int seqCompare = this.sequenceName.compareTo(o.sequenceName);
        if (seqCompare != 0) return seqCompare;
        int startCompare = Integer.compare(this.start, o.start);
        if (startCompare != 0) return startCompare;
        return Integer.compare(this.end, o.end);
    }

    public HashMap<String, Data> getData() {
        return data;
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

    /***
     * Order the samples in the data map according to the order in the header
     */
    public void alignSamplesWithHeader(String[] headerSamples) {
        HashMap<String, Data> alignedData = new LinkedHashMap<>();
        for (String sample : headerSamples) {
            alignedData.put(sample, this.data.getOrDefault(sample, null));
        }
        this.data = alignedData;
    }

    public String[] getSamples(){
        return data.keySet().toArray(new String[0]);
    }

    public double getScore(String sample){
        return data.get(sample).getScore();
    }

    public int getObservedKmers(String sample){
        return data.get(sample).getObservedKmers();
    }

    public int getVariations(String sample){
        return data.get(sample).getVariations();
    }

    public int getIbs(String sample){
        return data.get(sample).getIBS();
    }

    public void setIBS(String sample, int ibs){
        data.get(sample).setIBS(ibs);
    }

    public int getTotalKmers() {
        return totalKmers;
    }

    public void setWindowId(String windowId) {
        this.windowId = windowId;
    }

    public String toTSV(String sample) {
        return getWindowId() + "\t" + sequenceName + "\t" + start + "\t" + end + "\t" + getEffLength() + "\t" + totalKmers + "\t" + data.get(sample).toTSV();
    }

    public void setEffLength(int atgcCount) {
        this.effLength = atgcCount;
    }

    public int getEffLength() {
        return effLength;
    }

    public int getInnerDistance(String sample) {
        return data.get(sample).getInnerDistance();
    }

    public int getTailDistance(String sample) {
        return data.get(sample).getTailDistance();
    }

    public double getMeanKmerCount(String sample) {
        return data.get(sample).getMeanKmerCount();
    }
}
// EOF