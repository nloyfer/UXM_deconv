# UXM Reference Atlases

## Source

The reference atlases are derived from [Loyfer *et al.* (2022)](https://www.biorxiv.org/content/10.1101/2022.01.24.477547v1).
The `U25` atlases contain the top 25 specifically unmethylated marker blocks per cell type (~40 cell types, 900 markers total).
The `U250` atlases contain the top 250 markers per cell type.
Values represent the proportion of unmethylated reads (`U / (U + X)`) computed from WGBS pat files.

Original atlas files for **hg19** are as published. **hg38** and **chm13v2** versions were produced by liftover (see below).

---

## Liftover: hg19 → chm13v2.0\_maskedY\_rCRS

### Step 1 — Coordinate liftover

Genomic coordinates (`chr`, `start`, `end`) were lifted from hg19 to chm13v2.0\_maskedY\_rCRS using a genome coordinate liftover tool (e.g. CrossMap or UCSC liftOver with the appropriate chain file). Marker regions that could not be mapped are recorded in the corresponding `.unmapped.bed` files where present.

### Step 2 — Recalculate CpG indices

The `startCpG` and `endCpG` columns are genome-specific sequential CpG indices used by
[wgbstools](https://github.com/nloyfer/wgbs_tools) to look up reads in `.pat` files.
These indices are **not** transferable by coordinate liftover and must be recomputed for the target reference.

After liftover the correct indices were generated with:

```bash
# Write header
head -1 Atlas.U25.l4.chm13v2.tsv > Atlas.U25.l4.chm13v2.fixed.tsv

# Replace startCpG/endCpG with correct chm13v2 indices, then sort
# (wgbstools convert outputs: chr start end NEW_startCpG NEW_endCpG OLD_startCpG OLD_endCpG ...)
wgbstools convert \
    -L Atlas.U25.l4.chm13v2.tsv \
    --genome chm13v2.0_maskedY_rCRS \
    -p \
  | cut -f1-5,8- \
  | sort -k4,4n \
  >> Atlas.U25.l4.chm13v2.fixed.tsv

mv Atlas.U25.l4.chm13v2.fixed.tsv Atlas.U25.l4.chm13v2.tsv
```

The `sort -k4,4n` step is required because `wgbstools homog` expects markers sorted by `startCpG`
(monotonically increasing), whereas the liftover can change the relative ordering of markers across the genome.

### Notes

- The `name` column (format `chr:start-end`) retains coordinates from an intermediate liftover step
  and does not match the final `chr`/`start`/`end` columns. This is cosmetic only — `name` is used
  solely as a stable join key within UXM and its value does not affect results.
- 889 of the original 900 markers were successfully mapped to chm13v2; the remainder were lost during liftover.
