package nl.wur.bis.kcftools.Utils;

public enum Configs {
    KCF_VERSION("0.1"),
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
                    <ID=DI,Type=Integer,Description="Kmer Variation Distance">
                    <ID=SC,Type=Float,Description="Score">"""
    );

    private final String value;

    Configs(String value) {
        this.value = value;
    }

    public String getValue() {
        return value;
    }
}
