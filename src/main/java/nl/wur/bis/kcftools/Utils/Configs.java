package nl.wur.bis.kcftools.Utils;

import picocli.CommandLine;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;
/***
 * Enum to store the configuration values for the KCF file format
 */
public enum Configs {
    KCF_VERSION(getVersion()),
    KCF_SOURCE("kcftools"),
    KCF_DATE(HelperFunctions.getTodayDate()),
    KCF_INFO_LINES(
            """
                    <ID=EFFLEN,Type=Integer,Description="Effective length of the window">
                    <ID=IS,Type=Float,Description="Minimum score for the window">
                    <ID=XS,Type=Float,Description="Maximum score for the window">
                    <ID=MS,Type=Float,Description="Mean score for the window">
                    <ID=IO,Type=Integer,Description="Minimum observed kmers in the window">
                    <ID=XO,Type=Integer,Description="Maximum observed kmers in the window">
                    <ID=MO,Type=Integer,Description="Mean observed kmers in the window">
                    <ID=IV,Type=Integer,Description="Minimum variations in the window">
                    <ID=XV,Type=Integer,Description="Maximum variations in the window">
                    <ID=MV,Type=Integer,Description="Mean variations in the window">"""
    ),
    KCF_FORMAT_LINES(
            """
                    <ID=IB,Type=Integer,Description="IBS number">
                    <ID=VA,Type=Integer,Description="Variations">
                    <ID=OB,Type=Integer,Description="Observed kmers">
                    <ID=ID,Type=Integer,Description="Inner Distance">
                    <ID=LD,Type=Integer,Description="Kmer Variation Distance at the leftTail">
                    <ID=RD,Type=Integer,Description="Kmer Variation Distance at the rightTail">
                    <ID=KD,Type=Float,Description="Mean Kmer Depth">
                    <ID=SC,Type=Float,Description="Score">"""
    );

    private final String value;

    Configs(String value) {
        this.value = value;
    }

    public String getValue() {
        return value;
    }

    public static String getVersion() {
        try (InputStream in = Configs.class.getResourceAsStream("/version.properties")) {
            Properties props = new Properties();
            props.load(in);
            return props.getProperty("version", "UNKNOWN");
        } catch (IOException e) {
            return "UNKNOWN";
        }
    }

    public static class VersionProvider implements CommandLine.IVersionProvider {
        @Override
        public String[] getVersion() {
            return new String[]{ Configs.getVersion() };
        }
    }

}
//EOF