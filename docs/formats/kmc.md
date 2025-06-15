# KMC Database Format

This documentation is heavily inspired from the [KMC documentation](https://github.com/refresh-bio/KMC/blob/master/API.pdf)

KMC creates two files for each k-mer database:

- `.kmc_pre`: Contains prefix-related data and metadata.
- `.kmc_suf`: Contains suffixes and k-mer counters.

All integers are stored in **LSB (little-endian)** byte order.

---

## `.kmc_pre` File Structure

Contains the following sections:

1. `[marker]` — 4-byte string: `KMCP`
2. `[prefixes]`
3. `[map]`
4. `[header]`
5. `[header position]`
6. `[marker]` — End marker: `KMCP`

### Header Structure (`[header]`)

| Field                 | Type     | Description                                |
|----------------------|----------|--------------------------------------------|
| `kmer_length`        | `uint32` | Length of k-mers                           |
| `mode`               | `uint32` | 0 = count, 1 = quality-aware count         |
| `counter_size`       | `uint32` | Size of counters (1–4 bytes)               |
| `lut_prefix_length`  | `uint32` | Length of prefix removed from k-mers       |
| `signature_length`   | `uint32` | Length of k-mer signature                  |
| `min_count`          | `uint32` | Minimum count threshold                    |
| `max_count`          | `uint32` | Maximum count threshold                    |
| `total_kmers`        | `uint64` | Total number of k-mers stored              |
| `both_strands`       | `uchar`  | 1 = both strands used                      |
| `tmp[3]`             | `uchar[]`| Reserved                                   |
| `tmp[6]`             | `uint32[]`| Reserved                                  |
| `KMC_VER`            | `uint32` | KMC version (0x200 for KMC 2)              |

---

### Map Section (`[map]`)

- Array of `uint32` values.
- Size: `4^signature_length + 1`
- Used to locate the proper prefixes array in `[prefixes]`.

---

### Prefixes Section (`[prefixes]`)

- Contains multiple arrays of `uint64`, one for each signature.
- Each array has size `4^lut_prefix_length`.
- Followed by a guard `uint64`.

Used to locate k-mer suffix ranges in `.kmc_suf`.

---

## `.kmc_suf` File Structure

1. `[marker]` — 4-byte string: `KMCS`
2. `[data]` — Array of k-mer records
3. `[marker]` — End marker: `KMCS`

### Data Section (`[data]`)

- `record_t records[total_kmers]`
- Each record:
    - Suffix of the k-mer, packed (2 bits per base)
    - Counter value (1–4 bytes or 4-byte float)

---

!!! note
    - Use `[map]` and `[prefixes]` to find ranges in `.kmc_suf`.
    - Binary search can be applied within suffix ranges.
    - Only k-mers within min/max count thresholds are stored.

