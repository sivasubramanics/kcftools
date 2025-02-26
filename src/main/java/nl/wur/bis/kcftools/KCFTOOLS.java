package nl.wur.bis.kcftools;

import nl.wur.bis.kcftools.Plugins.*;
import nl.wur.bis.kcftools.Utils.Logger;
import picocli.CommandLine;
import picocli.CommandLine.*;

/***
 * This class is the main class for the KCFTOOLS tool
 */
@Command(name = "kcftools", description = "Tools to handle kmer counting based variations.",
        mixinStandardHelpOptions = true, version = "kcftools 0.1",
        subcommands = {
                GetVariants.class,
                Cohort.class,
                FindIBS.class,
                SplitKCF.class,
                GetAttributes.class,
                KCFToMatrix.class,
                KCFToTSV.class,
                CompareIBS.class
        })
public class KCFTOOLS {
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        final KCFTOOLS kcftools = new KCFTOOLS();
        final CommandLine cmd = new CommandLine(kcftools);
        Logger.info(KCFTOOLS.class.getSimpleName(), "Start processing");
        int exitCode = cmd.execute(args);
        long endTime = System.currentTimeMillis();
        String executionTime = String.format("%02d:%02d:%02d", (endTime - startTime) / 3600000,
                ((endTime - startTime) % 3600000) / 60000, ((endTime - startTime) % 60000) / 1000);
        switch (exitCode) {
            case 0:
                Logger.info(KCFTOOLS.class.getSimpleName(), "Total execution time: " + executionTime);
                break;
            case 2:
                break;
            default:
                Logger.error(KCFTOOLS.class.getSimpleName(), "Execution failed");
                break;
        }
        System.exit(exitCode);
    }
}

