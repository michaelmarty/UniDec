# Protein mixture MRS fixture

`protein_mix_two_injections.u8.zlib` contains 1,176 unsigned byte counts
in acquisition order, compressed with zlib level 9 to 537 bytes.
SHA-256 is `d0265f6a6a07d1b06cd15e11f1880cce5b4eac4c5f9c0693b796d782380ddea0`.
The data are unpublished. The user confirmed permission to redistribute
this derived fixture in UniDec.

The source is the lab paper archive under
`BHT/Temp/Paper/DMT/Protein Mix 21 Injections Set 2`. The first two RAW
and DMT pairs have the stem
`20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set2` and suffixes
`_1` and `_2`. The corresponding merged files and exports are under
`merged`. The acquisition name identifies a 0.1 mg/mL protein mixture on
an Agilent C8 column. The first RAW records positive-mode FTMS full scans
from m/z 500 to 3000. All 21 RAWs contain scans 1–588, and the first scan
and last scan of run 1 are at 0.0022250317 and 5.5025247 min. All merged
scan numbers and event times were checked against the RAW and DMT pairs.
For run number `r`, the merged scan is the local scan plus `588 * (r - 1)`,
and the merged time is the local RAW time plus the preceding RAW end times.

Count one ion per DMT row strictly inside m/z 673.6048583984375–2915.44970703125
and calibrated charge 1–36, using charge = intensity / 112670.893792227.
Count by the one-based merged scan number and fill empty scans with zero.
Keep scans 1–1,176 from the two-injection merge. Two complete cycles are
enough to exercise pooling while preserving the native bins and competing
peaks. The saved settings specify MRS, zero skipped cycles, and no scan
compression. The fixture uses no normalized intensities.

Direct addition of the two 588-scan count vectors reproduced the archived
two-injection MRS export and yields a unique apex at zero-based bin 385
(94 ions, 3.6098639 min in that export). The regression constructs a
linear time axis from the first RAW run's endpoints and allows 0.00017 min
for the largest measured axis deviation.

The independent retention-time reference is the separately acquired
ensemble TIC used in the unpublished Figure 2. Re-extracting the vendor
TIC from `20260327_BHT_agilentc8column_protienmix_ensemble_0o1mgmLs_1.raw`
in `BHT/Temp/Paper/Ensemble Protein Mix Triplicate` reproduced all 681
time/intensity pairs in `UniDec_Figures_and_Files/Protein mix Rep 1 chrom.txt`.
Its apex is RAW scan 445 at 3.6001669 min. The regression allows about
0.01746 min, the sum of one scan interval from each acquisition, to
account for sampling on the two time grids.
