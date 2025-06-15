# Limitations of the current implementation

This page outlines known limitations and important usage notes for `kcftools`, particularly when using the `getVariations` plugin. Please read these carefully to ensure correct and efficient usage of the tool.

---

## 1. KMC DB Compatibility

The `getVariations` plugin is compatible **only** with KMC databases generated using [`kmc`](https://github.com/refresh-bio/KMC) version **3.0.0 or newer**. Using databases created with older versions of KMC may result in unexpected behavior or incompatibility errors.

Additionally, support is currently **limited to KMC databases with a signature length of exactly 9**, specified using the `-p 9` parameter during database construction.

- **Files generated with different signature lengths are not supported.**
- Using incompatible databases may result in incorrect output, runtime errors, or program crashes.

It is **strongly recommended** to regenerate your KMC database using the appropriate version and settings to ensure compatibility with `kcftools`.

---

## 2. Memory Usage with `--memory` or `-m` Option

The `getVariations` plugin supports an optional `--memory` (or `-m`) flag that significantly improves runtime performance by loading the entire KMC database into system memory.

While this mode can dramatically accelerate computation, it **requires sufficient available RAM** and can lead to Java heap space errors on larger databases.

### To avoid memory-related crashes:

- **Increase the Java heap size manually** using the `-Xmx` JVM option:
    
        $ kcftools -Xmx16G getVariations ...

- **Alternatively**, define a default heap size by setting the environment variable `KCFTOOLS_HEAP_SIZE`:


        # Set the heap size to 16GB (or adjust as needed)
        $ export KCFTOOLS_HEAP_SIZE=16G

This will automatically configure the heap allocation when running `kcftools` with the recommended shell alias:

    $ alias kcftools='java -Xmx$KCFTOOLS_HEAP_SIZE -jar $PATH/kcftools.jar'


> Setting an appropriate heap size is crucial when working with large KMC databases or when using the `--memory` flag.

---

## Summary

| Limitation           | Description                                  |
| -------------------- | -------------------------------------------- |
| **KMC Version**      | Only `kmc` ≥ 3.0.0 supported                 |
| **Signature Length** | Only `-p 9` signature supported              |
| **Memory Mode**      | May cause heap errors without sufficient RAM |
| **Java Heap Size**   | Must be adjusted for large datasets          |

Please check this page regularly for updates regarding expanded support and future enhancements.

---

