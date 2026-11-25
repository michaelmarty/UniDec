# GDH FT-IM fixture

`gdh_32plus_counts.u8.zlib` holds one extracted ion count for each of the
166 frequency steps in a 4 kV bovine glutamate dehydrogenase (GDH) sweep.
The 93-byte file decompresses to 166 unsigned bytes. Its SHA-256 is
`9f1c3dc0caefbacfc829a79f2cefecd373e77cc11e7ce356b4075d0ad6418106`.

## Source and acquisition

The archive run has stem `03082024_GDH_4kV_stepsize3_repeate15_5to500_1` in
`CCS troubleshooting`
under `FT IM CD MS`. The vendor `.RAW` and the ion records in
`03082024_GDH_4kV_stepsize3_repeate15_5to500_1_2024-03-14-09-04-57.dmt`
are listed under `other/FT-IM-CD-MS/GDH/Varying Drift Voltage/` in
[MassIVE MSV000095092](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000095092)
(DOI `10.25345/C59W0996X`, CC0 1.0). The public DMT and archive DMT each
contain 37,872 ion rows and produce identical counts in every fixture bin.

The dual ion gates advance on the Orbitrap split-lens trigger. The run
uses 15 scans per step with gating frequencies 5, 8, …, 500 Hz, consistent
with the filename, saved `CDScanCompress 15` setting, and 2,490 scans in
the public RAW. RAW scan timestamps run from 0.014156 to 39.948498 min
with a median interval of 0.016044 min. Each interval differs from that
median by at most 0.0000021 min. The vendor header reports 800 ms DMT
injection and resolution 100,000 at *m/z* 400. The test uses a uniformly
spaced time axis of 0.240660 min per frequency step.
[Measurement and trigger method](https://doi.org/10.1021/acs.analchem.4c03379).

Scans 1–15 go into the first step. Empty steps retain a zero count.
`setup_ft` computes the sweep rate from the frequency span and the full
acquisition duration. The test passes an end-exclusive frequency of
503 Hz to represent 166 intervals of 3 Hz. The saved `FTend 500` is the
frequency of the last sampled step.

## Extraction and reference

Read `Ion` from the DMT SQLite database with `10484 <= Mz < 10584` around
the published 32+ signal at *m/z* 10534. Use the saved STORI charge
calibration of 53,000 slope units per charge to select
`31.5 * 53000 < Slope < 32.5 * 53000`. Count the 481 selected ions by
`(ScanNumber - 1) // 15`, with a maximum of 18 per bin. Write the 166
counts as `uint8` and compress them with `zlib` at level 9. The fixture
keeps the complete one-dimensional sweep and uses its uniform time grid
instead of storing time or frequency columns.

The 4 kV GDH 32+ trace in Figure 2 of Juetten et al. has a plotted apex
of approximately 128.5 nm². The publication identifies the condition
and charge state but gives no filename for the representative trace.
Using Mason–Schamp with the archived configuration (0.104 m between
gates, 2,000 V across them, 25 °C, 760 torr, nitrogen at 28.0134 Da, and
zero time offset), the plotted apex corresponds to about 98.23 ms.
[Peer-reviewed reference](https://doi.org/10.1021/acs.analchem.4c03379) ·
[Author manuscript with Figure 2](https://doi.org/10.26434/chemrxiv-2024-cbqj1).

An independent calculation from the fixture counts used a half-Hann
window, 15-point local cubic least-squares detrending, and direct sums of
the 256-point Fourier harmonics. Its maximum was bin 76 (zero-based) at
98.958333 ms, which the engine reproduced. The test requires that bin
and allows one native bin (`1000 / (166 * 3)`, about 2.008 ms) when
comparing with the plotted peak. This is the resolution of the acquired
sweep before zero padding and covers the precision of reading the plot.
