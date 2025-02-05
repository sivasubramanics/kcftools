package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.HelperFunctions;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Objects;

/***
 * Class to represent a window in the KCF file
 */
public class Window implements Comparable<Window> {
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
    float minScore;
    float maxScore;
    float meanScore;


    public Window(String windowId, String sequenceName, int start, int end){
        this.windowId = windowId;
        this.sequenceName = sequenceName;
        this.start = start;
        this.end = end;
        this.effLength = 0;
        this.totalKmers = 0;
        this.data = new LinkedHashMap<>();
    }

    public Window(String[] fields, String[] samples){
        this.windowId = fields[3];
        this.sequenceName = fields[0];
        this.start = Integer.parseInt(fields[1]);
        this.end = Integer.parseInt(fields[2]);
        this.totalKmers = Integer.parseInt(fields[4]);
        this.effLength = Integer.parseInt(getInfoFieldMap(fields[5]).get("EFFLEN"));
        this.data = new LinkedHashMap<>();
        for (int i = 7; i < fields.length; i++){
            String[] sampleData = fields[i].split(":");
            // sampleData: ibs:variations:observedKmers:innerDistance:tailDistance:score
            data.put(samples[i-7], new Data(Integer.parseInt(sampleData[2]), Integer.parseInt(sampleData[1]), Integer.parseInt(sampleData[3]), Integer.parseInt(sampleData[4]), totalKmers, effLength, sampleData[0]));
        }
    }

    public void addData(HashMap<String, Data> data){
        for (String sample : data.keySet()){
            if (this.data.containsKey(sample)){
                HelperFunctions.log("error", "Window", "Sample " + sample + " already exists in window " + windowId);
            } else {
                this.data.put(sample, data.get(sample));
            }
        }
    }

    public synchronized void addData(String sample, int observedKmers, int variations, int innerDistance, int tailDistance, String ibs) {
        Data d = data.computeIfAbsent(sample, k -> new Data(0, 0, 0, 0, totalKmers, effLength));
        d.observedKmers = observedKmers;
        d.variations = variations;
        d.innerDistance = innerDistance;
        d.tailDistance = tailDistance;
        d.score = d.computeScore(totalKmers, effLength);
        d.ibs = "N".equals(ibs) ? -1 : Integer.parseInt(ibs);
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
        return "GT:VA:OB:ID:TD:SC";
    }

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

    public float getScore(String sample){
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
        // return sequenceName + "\t" + start + "\t" + end + "\t" + totalKmers + "\t" + data.get(sample).toTSV();
        return sequenceName + "\t" + start + "\t" + end + "\t" + totalKmers + "\t" + data.get(sample).toTSV();
    }

    public void setEffLength(int atgcCount) {
        this.effLength = atgcCount;
    }

    public int getEffLength() {
        return effLength;
    }
}
// EOF