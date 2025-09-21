package nl.wur.bis.kcftools.Data;

import nl.wur.bis.kcftools.Utils.Configs;
import nl.wur.bis.kcftools.Utils.HelperFunctions;
import nl.wur.bis.kcftools.Utils.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * KCFHeader class to store the header information of the KCF file (still there is hope to improve and make it more flexible)
 */
public class KCFHeader implements Comparable<KCFHeader> {

    private final String CLASSNAME = this.getClass().getSimpleName();
    private String version;
    private String source;
    private final String date;
    private String reference;
    private LinkedHashMap<String, Integer> contigs;
    private final String infoLines;
    private final String formatLines;
    private List<String> commandLines;
    private String[] samples;
    // window [0], step [1], kmer [2], IBS [3], numWindow [4], weightInnerDistance [5], weightTailDistance [6], weightKmerRatio [7] as parameters in the order
    private final Pair[] params = new Pair[8];

    public KCFHeader(){
        this.version = Configs.KCF_VERSION.getValue();
        this.source = Configs.KCF_SOURCE.getValue();
        this.date = HelperFunctions.getTodayDate();
        this.reference = "";
        this.contigs = null;
        this.infoLines = Configs.KCF_INFO_LINES.getValue();
        this.formatLines = Configs.KCF_FORMAT_LINES.getValue();
    }

    public KCFHeader(String reference, LinkedHashMap<String, Integer> contigs){
        this.version = Configs.KCF_VERSION.getValue();
        this.source = Configs.KCF_SOURCE.getValue();
        this.date = HelperFunctions.getTodayDate();
        this.reference = reference;
        this.contigs = contigs;
        this.infoLines = Configs.KCF_INFO_LINES.getValue();
        this.formatLines = Configs.KCF_FORMAT_LINES.getValue();
    }

    public KCFHeader(String headerLines){
        this.version = Configs.KCF_VERSION.getValue();
        this.source = Configs.KCF_SOURCE.getValue();
        this.date = HelperFunctions.getTodayDate();
        this.reference = "";
        this.contigs = null;
        this.infoLines = Configs.KCF_INFO_LINES.getValue();
        this.formatLines = Configs.KCF_FORMAT_LINES.getValue();
        String[] lines = headerLines.split("\n");
        for (String line : lines){
            if (line.startsWith("##reference=")){
                this.reference = line.substring(12);
            } else if (line.startsWith("##contig=")){
                String[] contigLine = line.substring(10, line.length()-1).split(",");
                String contigName = contigLine[0].substring(3);
                int contigLength = Integer.parseInt(contigLine[1].substring(7));
                addContig(contigName, contigLength);
            } else if (line.startsWith("##CMD=")){
                addCommandLine(line.substring(6));
            } else if (line.startsWith("#CHROM")){
                // split line by tab and get sample names from the 6th index
                String[] fields = line.split("\t");
                samples = new String[fields.length - 7];
                System.arraycopy(fields, 7, samples, 0, fields.length - 7);
            } else if (line.startsWith("##PARAM=")){
                Pair param = getParam(line);
                switch (param.getKey()){
                    case "window":
                        params[0] = param;
                        break;
                    case "step":
                        params[1] = param;
                        break;
                    case "kmer":
                        params[2] = param;
                        break;
                    case "IBS":
                        params[3] = param;
                        break;
                    case "nwindow":
                        params[4] = param;
                        break;
                    case "wti":
                        params[5] = param;
                        break;
                    case "wtt":
                        params[6] = param;
                        break;
                    case "wtk":
                        params[7] = param;
                        break;
                }
            }
        }
    }

    /***
     * Get a String array of all the contig/chromosome/scaffold names
     * @return
     */
    public String[] getContigs(){
        return contigs != null ? contigs.keySet().toArray(new String[0]) : null;
    }

    /***
     * Get contigs ID (just an integer index) for the contig name
     */
    public int getContigID(String contigName){
        if (contigs == null || !contigs.containsKey(contigName)){
            Logger.error(CLASSNAME, "Contig " + contigName + " not found in the KCF header");
            return -1;
        }
        return new ArrayList<>(contigs.keySet()).indexOf(contigName);
    }

    /***
     * Get a set of Param Key-Value pairs
     */
    private Pair getParam(String line) {
        String[] fields = line.substring(9, line.length()-1).split(",");
        String key = fields[0].substring(3);
        String value = fields[1].substring(6);
        return new Pair(key, value);
    }

    /***
     * Get Weight for obs/total kmers
     */
    public double getWeightKmerRatio(){
        return this.params[7] != null ? Double.parseDouble(this.params[7].getValue()) : 0.0;
    }

    /***
     * Set Weight for obs/total kmers
     */
    public void setWeightKmerRatio(double wtk){
        this.params[7] = new Pair("wtk", String.valueOf(wtk));
    }

    /***
     * Get Weight for tail distance
     */
    public double getWeightTailDist(){
        return this.params[6] != null ? Double.parseDouble(this.params[6].getValue()) : 0.0;
    }

    /***
     * Set Weight for tail distance
     */
    public void setWeightTailDist(double wtt){
        this.params[6] = new Pair("wtt", String.valueOf(wtt));
    }

    /***
     * Get Weight for inner distance
     */
    public double getWeightInnerDist(){
        return this.params[5] != null ? Double.parseDouble(this.params[5].getValue()) : 0.0;
    }

    /***
     * Set Weight for inner distance
     */
    public void setWeightInnerDist(double wti){
        this.params[5] = new Pair("wti", String.valueOf(wti));
    }

    /***
     * Get the number of windows
     */
    public int getWindowCount(){
        return this.params[4] != null ? Integer.parseInt(this.params[4].getValue()) : 0;
    }

    /***
     * Set the number of windows (this is important for the IBS processing)
     */
    public void setWindowCount(int nunWindow){
        this.params[4] = new Pair("nwindow", String.valueOf(nunWindow));
    }

    public boolean isIBS(){
        return this.params[3] != null && Boolean.parseBoolean(this.params[3].getValue());
    }

    public void setIBS(boolean isIBS){
        this.params[3] = new Pair("IBS", String.valueOf(isIBS));
    }

    public int getWindowSize(){
        return this.params[0] != null ? Integer.parseInt(this.params[0].getValue()) : 0;
    }

    public int getStepSize(){
        return this.params[1] != null ? Integer.parseInt(this.params[1].getValue()) : 0;
    }

    public int getKmerSize(){
        return this.params[2] != null ? Integer.parseInt(this.params[2].getValue()) : 0;
    }

    public void setWindowSize(int windowSize){
        this.params[0] = new Pair("window", String.valueOf(windowSize));
    }

    public void setStepSize(int stepSize) {
        this.params[1] = new Pair("step", String.valueOf(stepSize));
    }

    public void setKmerSize(int kmerSize){
        this.params[2] = new Pair("kmer", String.valueOf(kmerSize));
    }

    public String[] getSamples(){
        return samples;
    }

    public void setSample(String[] sampleNames){
        this.samples = sampleNames;
    }

    /***
     * Add a sample name to the list of samples
     */
    public void addSample(String sampleName){
        if (samples == null){
            samples = new String[0];
        }
        String[] newSampleNames = new String[samples.length + 1];
        System.arraycopy(samples, 0, newSampleNames, 0, samples.length);
        newSampleNames[samples.length] = sampleName;
        samples = newSampleNames;
    }

    /***
     * Append a list of sample names to the list of samples
     */
    public void addSample(String[] sampleNames){
        if (samples == null){
            samples = new String[0];
        }
        String[] newSampleNames = new String[samples.length + sampleNames.length];
        System.arraycopy(samples, 0, newSampleNames, 0, samples.length);
        System.arraycopy(sampleNames, 0, newSampleNames, samples.length, sampleNames.length);
        samples = newSampleNames;
    }

    /***
     * Add a command line to the list of command lines
     */
    public void addCommandLine(String commandLine){
        if (commandLines == null){
            commandLines = new ArrayList<>();
        }
        commandLines.add(commandLine);
    }

    /***
     * Add a contig name and length to the list of contigs
     */
    public void addContig(String name, int length){
        if (contigs == null){
            contigs = new LinkedHashMap<>();
        }
        contigs.put(name, length);
    }

    /***
     * Reference genome of the which the getVariants are called
     */
    public void setReference(String reference){
        this.reference = reference;
    }

    public void writeHeader(BufferedWriter writer) throws IOException {
        writer.write(toString());
    }

    public void setVersion(String v) {
        this.version = v;
    }

    public void setSource(String s) {
        this.source = s;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("##format=KCF").append(version).append("\n");
        sb.append("##date=").append(date).append("\n");
        sb.append("##source=").append(source).append("\n");
        sb.append("##reference=").append(reference).append("\n");
        if (contigs != null) {
            for (String contig : contigs.keySet()) {
                sb.append("##contig=<ID=").append(contig).append(",length=").append(contigs.get(contig)).append(">\n");
            }
        }
        if (infoLines != null) {
            for (String infoLine : infoLines.split("\n")) {
                sb.append("##INFO=").append(infoLine).append("\n");
            }
        }
        if (formatLines != null) {
            for (String formatLine : formatLines.split("\n")) {
                sb.append("##FORMAT=").append(formatLine).append("\n");
            }
        }
        for (Pair param : params){
            if (param != null){
                sb.append("##PARAM=<ID=").append(param.getKey()).append(",value=").append(param.getValue()).append(">\n");
            }
        }
        if (commandLines != null) {
            for (String commandLine : commandLines) {
                sb.append("##CMD=").append(commandLine).append("\n");
            }
        }
        sb.append("#CHROM\tSTART\tEND\tID\tTOTAL_KMERS\tINFO\tFORMAT");
        if (samples != null){
            for (String sample : samples){
                sb.append("\t").append(sample);
            }
        }
        sb.append("\n");
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        KCFHeader kcfHeader = (KCFHeader) o;
        if (getWindowSize() != kcfHeader.getWindowSize()) {
            Logger.error(CLASSNAME, "Window size mismatch between the KCFs");
            return false;
        }
        if (getKmerSize() != kcfHeader.getKmerSize()) {
            Logger.error(CLASSNAME, "Kmer size mismatch between the KCFs");
            return false;
        }
        if (isIBS() != kcfHeader.isIBS()) {
            Logger.error(CLASSNAME, "IBS processing mismatch between the KCFs");
            return false;
        }
        if (getWindowCount() != kcfHeader.getWindowCount()) {
            Logger.error(CLASSNAME, "Number of windows mismatch between the KCFs");
            return false;
        }
        if (getWeightInnerDist() != kcfHeader.getWeightInnerDist()) {
            Logger.error(CLASSNAME, "Weight Inner Distance mismatch between the KCFs");
            return false;
        }
        if (getWeightTailDist() != kcfHeader.getWeightTailDist()) {
            Logger.error(CLASSNAME, "Weight Tail Distance mismatch between the KCFs");
            return false;
        }
        if (getWeightKmerRatio() != kcfHeader.getWeightKmerRatio()) {
            Logger.error(CLASSNAME, "Weight Kmer Ratio mismatch between the KCFs");
            return false;
        }
        if (getStepSize() != kcfHeader.getStepSize()) {
            Logger.error(CLASSNAME, "Step size mismatch between the KCFs");
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        return Objects.hash(getWindowSize(), getKmerSize(), isIBS(), getWindowCount());
    }

    @Override
    public int compareTo(@NotNull KCFHeader o) {
        int winCompare = Integer.compare(this.getWindowSize(), o.getWindowSize());
        if (winCompare != 0) {
            Logger.error(CLASSNAME, "Window size mismatch between the KCFs");
            return winCompare;
        }
        int kmerCompare = Integer.compare(this.getKmerSize(), o.getKmerSize());
        if (kmerCompare != 0) {
            Logger.error(CLASSNAME, "Kmer size mismatch between the KCFs");
            return kmerCompare;
        }
        int nwinCompare = Integer.compare(this.getWindowCount(), o.getWindowCount());
        if (nwinCompare != 0) {
            Logger.error(CLASSNAME, "Number of windows mismatch between the KCFs");
            return nwinCompare;
        }
        int ibsCompare = Boolean.compare(this.isIBS(), o.isIBS());
        if (ibsCompare != 0) {
            Logger.error(CLASSNAME, "IBS processing mismatch between the KCFs");
            return ibsCompare;
        }
        int wtiCompare = Double.compare(this.getWeightInnerDist(), o.getWeightInnerDist());
        if (wtiCompare != 0) {
            Logger.error(CLASSNAME, "Weight Inner Distance mismatch between the KCFs");
            return wtiCompare;
        }
        int wttCompare = Double.compare(this.getWeightTailDist(), o.getWeightTailDist());
        if (wttCompare != 0) {
            Logger.error(CLASSNAME, "Weight Tail Distance mismatch between the KCFs");
            return wttCompare;
        }
        int wtkCompare = Double.compare(this.getWeightKmerRatio(), o.getWeightKmerRatio());
        if (wtkCompare != 0) {
            Logger.error(CLASSNAME, "Weight Kmer Ratio mismatch between the KCFs");
            return wtkCompare;
        }
        return 0;
    }

    /***
     * Merge the header of the KCF file with the current header
     */
    public void mergeHeader(KCFHeader tmpHeader) {
        if (!this.equals(tmpHeader)){
            Logger.error(CLASSNAME, "Headers mismatch found in the KCF files");
        }
        if (tmpHeader.getSamples() != null){
            addSample(tmpHeader.getSamples());
        }
        if (tmpHeader.getCMDlines() != null){
            for (String cmd : tmpHeader.getCMDlines()){
                addCommandLine(cmd);
            }
        }
    }

    /***
     * Get the command lines used to generate the KCF file
     */
    private String[] getCMDlines() {
        return commandLines != null ? commandLines.toArray(new String[0]) : null;
    }

    /***
     * Check if the sample name is in the list of samples
     */
    public boolean hasSample(String sampleName) {
        return samples != null && Arrays.asList(samples).contains(sampleName);
    }

    /***
     * Get the weight of the parameters
     */
    public double[] getWeights() {
        return new double[]{getWeightInnerDist(), getWeightTailDist(), getWeightKmerRatio()};
    }
}

/**
 * Pair class to store key-value pairs
 */
class Pair{
    private final String key;
    private final String value;

    public Pair(String key, String value){
        this.key = key;
        this.value = value;
    }

    public String getKey(){
        return key;
    }

    public String getValue(){
        return value;
    }
}
//EOF