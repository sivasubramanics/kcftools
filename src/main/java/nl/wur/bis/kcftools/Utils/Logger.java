package nl.wur.bis.kcftools.Utils;

import java.text.SimpleDateFormat;
import java.util.Date;

public class Logger {
    public enum LogLevel {
        INFO, WARNING, ERROR, DEBUG, UNKNOWN
    }

    public static void log(LogLevel logLevel, String className, String message) {
        // assuming max className is 20 characters long
        if (className.length() < 20) {
            className = String.format("%-20s", className);
        }

        String timeString = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss:SSS").format(new Date());
        String logPrefix = timeString + " - " + logLevelToString(logLevel) + " - " + className + " - ";

        if (logLevel == LogLevel.DEBUG && System.getProperty("DEBUG") == null) {
            return;
        }

        System.out.println(logPrefix + message);

        if (logLevel == LogLevel.ERROR) {
            System.exit(1);
        }
    }

    private static String logLevelToString(LogLevel logLevel) {
        return switch (logLevel) {
            case INFO -> "INFO    ";
            case WARNING -> "WARNING ";
            case ERROR -> "ERROR   ";
            case DEBUG -> "DEBUG   ";
            default -> "UNKNOWN ";
        };
    }

    public static void info(String className, String message) {
        log(LogLevel.INFO, className, message);
    }

    public static void warning(String className, String message) {
        log(LogLevel.WARNING, className, message);
    }

    public static void error(String className, String message) {
        log(LogLevel.ERROR, className, message);
    }

    public static void debug(String className, String message) {
        log(LogLevel.DEBUG, className, message);
    }
}
