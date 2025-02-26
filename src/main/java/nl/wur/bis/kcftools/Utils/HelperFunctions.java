package nl.wur.bis.kcftools.Utils;


import picocli.CommandLine;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;


/**
 * Utility Helper functions for Meganome
 */
public class HelperFunctions {

    private static final String[] DEPENDENCIES = {"kmc"};

    private static final String CLASS_NAME = HelperFunctions.class.getSimpleName();

    /**
     * decompresses a gzip file
     */
    public static String decompressGZipFile(String inFile) throws IOException {
        String fileName = inFile.replace(".gz", "");
        FileInputStream fis = new FileInputStream(inFile);
        GZIPInputStream gis = new GZIPInputStream(fis);
        FileOutputStream fos = new FileOutputStream(fileName);
        byte[] buffer = new byte[1024];
        int len;
        while((len = gis.read(buffer)) != -1){
            fos.write(buffer, 0, len);
        }
        gis.close();
        fos.close();
        return fileName;
    }

    /***
     * Execute a command and return the output
     */
    public static Process initProcess(String[] execCommand,
                                      String outputFile,
                                      String errorFile,
                                      String inputFile,
                                      Integer waitSeconds) {

        Process returnVal = null;
        String execCommandProcName = execCommand[0];

        ProcessBuilder builder = new ProcessBuilder(execCommand);
        if (outputFile != null) builder.redirectOutput(new File(outputFile));
        if (errorFile != null) builder.redirectError(new File(errorFile));
        if (inputFile != null) builder.redirectInput(new File(inputFile));

        try {
            returnVal = builder.start();
            if (waitSeconds != null) {
                returnVal.waitFor(waitSeconds, TimeUnit.SECONDS);
            } else {
                returnVal.waitFor();
            }
        } catch (Exception e) {
            System.out.println(execCommandProcName + "Exception in process" + e);
            System.out.println(execCommandProcName + "Error File Contents" + errorFile);
        }

        return returnVal;
    }

    /***
     * Execute a command and return the output
     */
    public static String[] makeExecString(String execString) {
        return execString.split(" ");
    }

    /***
     * Execute a command and return the output
     */
    public static String makeProcName(String[] execString) {

        String returnVal = execString[0];

        if (returnVal.equals("python")) {
            returnVal = execString[1];
        }
        return returnVal;
    }

    /***
     * Execute a command and return the output
     */
    public static boolean tryExec(String execString, String outputFile, String errorFile, String inputFile) {

        String[] execArray = makeExecString(execString);
        String executedProcName = makeProcName(execArray);

        Process p = initProcess(execArray, outputFile, errorFile, inputFile, null);

        if (p == null) {
            System.out.println(executedProcName + " Process did not execute correctly - no data returned");
            System.exit(1);
            return false;
        }
        if (p.exitValue() != 0) {
            String message = executedProcName + " exited with a return code of " + p.exitValue();
            System.out.println(executedProcName + " " + message + " " + errorFile);
            System.exit(1);
            return false;
        }

        return true;
    }

    /***
     * Execute a command and return the output
     */
    public static boolean tryExec(String toExec) {
        return tryExec(toExec, null, null);
    }

    /***
     * Execute a command and return the output
     */
    public static boolean tryExec(String execString, String outputFile, String errorFile) {
        return tryExec(execString, outputFile, errorFile, null);
    }

    /***
     * Execute a command and return the output
     */
    public static void checkExecutables(String[] executables) {
        for (String executable : executables) {
            if (!tryExec("which " + executable)) {
                throw new RuntimeException("Executable " + executable + " not found in path");
            }
        }
    }

    /***
     * Check the dependencies are installed for the meganome pipeline
     */
    public static void checkDependencies() {
        try{
            checkExecutables(DEPENDENCIES);
        }
        catch (Exception e){
            System.out.println(e.getMessage());
            System.out.println("All the dependencies are not installed. Please install them before running meganome. { "+ String.join(", ", DEPENDENCIES) + " }");
            System.exit(1);
        }
    }

    /***
     * Check if a directory exists
     */
    public static boolean checkDirectoryExists(String baseDir) {
        return checkFileExists(baseDir);
    }

    /***
     * Check if a file exists
     */
    public static boolean checkFileExists(String inFile) {
        return Files.exists(Path.of(inFile));
    }

    /***
     * Reverse complement a sequence
     */
    public static String reverseComplement(String sequence) {
        char[] sequenceArray = sequence.toCharArray();
        char[] reverseSequence = reverseComplement(sequenceArray);
        return new String(reverseSequence);
    }

    /***
     * Reverse complement a sequence
     */
    public static char[] reverseComplement(char[] sequence) {
        int length = sequence.length;
        char[] reverseSequence = new char[length];

        for (int i = 0; i < length; i++) {
            reverseSequence[length - i - 1] = getComplement(sequence[i]);
        }
        return reverseSequence;
    }

    /***
     * Get the complement of a nucleotide
     */
    private static char getComplement(char nucleotide) {
        return switch (nucleotide) {
            case 'A' -> 'T';
            case 'a' -> 't';
            case 'T' -> 'A';
            case 't' -> 'a';
            case 'C' -> 'G';
            case 'c' -> 'g';
            case 'G' -> 'C';
            case 'g' -> 'c';
            case 'N' -> 'N';
            case 'n' -> 'n';
            default -> throw new IllegalArgumentException("Invalid nucleotide: " + nucleotide);
        };
    }


    /***
     * Check if a tool is installed
     */
    public static boolean isInstalled(String toolName) {
        try {
            Process process = Runtime.getRuntime().exec("which " + toolName);
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.contains(toolName)) {
                    return true;
                }
            }
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error checking if " + toolName + " is installed");
        }
        return false;
    }

    public static boolean isCompressed(String fileName){
        return isCompressed(new File(fileName));
    }

    /***
     * Check if a file is compressed
     */
    public static boolean isCompressed(File fileName){
        try{
            InputStream fileStream = new FileInputStream(fileName);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            gzipStream.close();
        } catch(ZipException e){
            return false;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return true;
    }

    /***
     * Fold a sequence for specified length
     */
    public static String fold_seq(String seq, int len) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < seq.length(); i += len) {
            sb.append(seq, i, Math.min(i + len, seq.length()));
            sb.append("\n");
        }
        return sb.toString();
    }

    /***
     * Create a directory
     */
    public static void createDirectory(String inDir) { tryExec("mkdir -p " + inDir); }

    /***
     * Delete a directory
     * @param inDir: the directory to delete
     */
    public static void deleteDirectory(String inDir) {
        if (inDir.equals("/")) {
            throw new RuntimeException("Cannot delete root directory");
        }
        if (checkDirectoryExists(inDir)) {
            tryExec("rm -rf " + inDir);
        }
    }

    /***
     * Get the base name of a file
     */
    public static String getBaseName(String inFile) {
        return Path.of(inFile).getFileName().toString();
    }

    /***
     * Based on the total and current number show the print the progress on STDOUT, self deleting the last line
     */
    public static void showProgress(int totalNum, int currentNum) {
        String progress = String.format("%.2f", (currentNum * 100.0) / totalNum);
        System.out.print("\r" + "Processed: " + progress + "%");
        if (currentNum == totalNum) {
            System.out.print("\r");
        }
    }

    /***
     * This method is used to delete a file
     */
    public static void deleteFile(String inFile) {
        try {
            Files.delete(Path.of(inFile));
        } catch (IOException e) {
            Logger.error(CLASS_NAME, "Error deleting file: " + inFile);
        }
    }


    /***
     * This method is used to check if a kmer contains non ACGT characters
     */
    public static boolean containsNonACGT(String kmer) {
        char[] kmerChars = kmer.toCharArray();
        return containsNonACGT(kmerChars);
    }

    /***
     * This method is used to check if a kmer contains non ACGT characters
     */
    public static boolean containsNonACGT(char[] kmer) {
        for (char base : kmer) {
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                return true;
            }
        }
        return false;
    }

    /***
     * This method is used to compare two byte arrays
     */
    public static int compareByteArray(byte[] suffixArray, byte[] suffix) {
        int i;
        for(i = 0; i < suffixArray.length; ++i)
            if(suffixArray[i] != suffix[i])
                break;
        if(i == suffixArray.length)
            return 0;
        else if((suffixArray[i] & 0xFF) < (suffix[i] & 0xFF) )
            return -1;
        else
            return 1;
    }

    /***
     * This method is to compare two files, and return true if first file is older than the second file (last modifed time)
     */
    public static boolean isOlder(File fileOne, File fileTwo) {
        return fileOne.lastModified() < fileTwo.lastModified();
    }

    /***
     * This method returns available RAM size based on the unity provided. If the unit is gb then return in gb, if mb then return in mb
     */
    public static String getMaxRAMSize(String unit) {
        // maximum avaiable memory in bytes
        long totalMemoryBytes = Runtime.getRuntime().maxMemory();
        long availableMemoryBytes = Runtime.getRuntime().freeMemory();
        double totalMemory;
        double availableMemory;

        if (unit == null) {
            return String.format("%d:%d", totalMemoryBytes, availableMemoryBytes);
        }

        switch (unit.toUpperCase()) {
            case "GB":
                totalMemory = totalMemoryBytes / (1024.0 * 1024.0 * 1024.0);
                availableMemory = availableMemoryBytes / (1024.0 * 1024.0 * 1024.0);
                break;
            case "MB":
                totalMemory = totalMemoryBytes / (1024.0 * 1024.0);
                availableMemory = availableMemoryBytes / (1024.0 * 1024.0);
                break;
            case "KB":
                totalMemory = totalMemoryBytes / 1024.0;
                availableMemory = availableMemoryBytes / 1024.0;
                break;
            default:
                System.out.println("Invalid unit. Please use GB or MB. Returning in bytes.");
                totalMemory = totalMemoryBytes;
                availableMemory = availableMemoryBytes;
        }
        return String.format("%.2f:%.2f", totalMemory, availableMemory);
    }

    /***
     * This method is used to check if a string is numeric
     */
    public static boolean isNumeric(String string) {
        try {
            Integer.parseInt(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    /***
     * This method is used to get the current date
     */
    public static String getTodayDate() {
        return new SimpleDateFormat("yyyy-MM-dd").format(new Date());
    }

    /***
     * This method is used to get the command line used to run the program
     */
    public static String getCommandLine() {
        return String.join(" ", System.getProperty("sun.java.command"));
    }

    /***
     * Print the command line options
     */
    public static void printCommandLine(CommandLine commandLine, String className) {
        String name;
        String value;
        Logger.info(className, "========== CMD options - " + className + " ==========");
        for (CommandLine.Model.OptionSpec option : commandLine.getCommandSpec().options()) {
            if (option.isOption()) {
                // get long name of the command
                if (option.names().length == 1) {
                    name = option.names()[0];
                }
                else {
                    name = option.names()[1];
                }
                // if getValue is null skip
                if (option.getValue() == null) {
                    continue;
                }
                value = option.getValue().toString();
                Logger.info(className, String.format("%-15s: %s", name, value));
            }
        }
        Logger.info(className, "==================================================");
    }

    /***
     * Print the memory usage
     */
    public static void printMaxMemoryUsage() {
        Runtime runtime = Runtime.getRuntime();
        long maxMemory = runtime.maxMemory();
        long allocatedMemory = runtime.totalMemory();
        long freeMemory = runtime.freeMemory();
        long usedMemory = allocatedMemory - freeMemory;

        Logger.info(CLASS_NAME, "============= Memory Usage Statistics ============");
        Logger.info(CLASS_NAME, String.format("%-25s: %.2f", "Max Memory (GB)", maxMemory / (1024.0 * 1024 * 1024)));
        Logger.info(CLASS_NAME, String.format("%-25s: %.2f", "Allocated Memory (GB)", allocatedMemory / (1024.0 * 1024 * 1024)));
        Logger.info(CLASS_NAME, String.format("%-25s: %.2f", "Free Memory (GB)", freeMemory / (1024.0 * 1024 * 1024)));
        Logger.info(CLASS_NAME, String.format("%-25s: %.2f", "Used Memory (GB)", usedMemory / (1024.0 * 1024 * 1024)));
        Logger.info(CLASS_NAME, "==================================================");

    }
}
// EOF