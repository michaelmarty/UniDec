# β-galactosidase HT-SEC TIC

`bgal_5bit_tic.u8` is the total ion chromatogram (TIC) for the first
5-bit β-galactosidase-only run in Figure 2 of Sanders et al. It contains
474 unsigned-byte ion counts in scan-bin order, with no header. The grid
is implicit. For zero-based array index `i`, time is `(i + 1) * 41 / 474`
minutes. SHA-256 is
`b5f2236a1d42c39ed167992138fd2be568bf7d05ee8f1b4ddcf32e6ca878eb73`.
The paper identifies MassIVE MSV000095091 as the source of its raw data.
MassIVE lists the dataset under CC0 1.0. The fixture contains only ion
counts derived from that experimental run.
[Published study](https://doi.org/10.1021/acs.analchem.4c03248) ·
[Methods and Figure 2](https://doi.org/10.26434/chemrxiv-2024-h08x5-v2) ·
[MassIVE dataset and license](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=87c596b5bd0646e0ac5237699b774ba4).

## Extraction and acquisition mapping

The paper archive contains `Bgal bit len/5 bit/20240411 Bgal bit5 zp10 1.RAW.zip`
and the companion `20240411 Bgal bit5 zp10 1_2024-04-26-03-44-37_unidecfiles.zip`.
Use the `data` array from the companion's `*_rawdata.npz`. Keep events with
`5000.62109375 < m/z < 17969.73828125` and STORI slopes strictly between
`25 * 116000` and `56 * 116000`, using the saved charge filter. Count events
by `scan // 10`, then write bins 1 through 474 as unsigned bytes. All 474
counts match the archived `TIC_chrom.txt` exactly. Including the four ions
in bin 0 gives 22,725 ions, the count for 5-bit replicate 1 in Table 1.

The Q Exactive HF UHMR ran in DMT mode with 200 ms ion injection.
The LC valve injected for 1.5 s at the start of each occupied one-minute
slot. The controller triggered MS acquisition at sequence start. The
31-slot order is `0000100101100111110001101110101`, followed by ten
non-injection slots. The saved kernel agrees with this order and
the saved processing settings (`htbit=5`, `HTtimepad=10`, `HTtimeshift=10`,
`HTcycleindex=-1`, `HTksmooth=1`, scan compression 10). Independently of
that kernel, summing copies of control 1 shifted by the programmed
injection times gives a Pearson correlation of 0.973 with the raw TIC
at the documented origin, versus 0.512 and 0.537 at offsets of one
cycle in either direction (using a circular shift of 12 bins).

Thermo RAW timestamps span 0.003612 to 41.003125 min across 4,743 scans,
with a median interval of 0.008646 min. Every scan assigned to bins 1
through 474 lies within one compressed bin (0.086498 min) of its
archived grid time. The regression uses this archived uniform grid.

## Independent reference and tolerance

Use the separate 12-minute, single-injection traces shown in Figure 2
of the peer-reviewed paper. The corresponding `TIC_chrom.txt` files in
the archive's `Bgal bit len/1 bit` analysis ZIPs give these apex times.
The final column gives the first and last vendor RAW timestamp assigned
to each apex bin.

| Run | Archived apex (min) | RAW apex-bin span (min) |
| --- | ---: | ---: |
| `20240411 Bgal bit1 1` | 9.684211 | 9.678530–9.756344 |
| `20240411 Bgal bit1 2` | 9.597744 | 9.592070–9.669885 |
| `20240411 Bgal bit1 3` | 9.597744 | 9.592070–9.669884 |

The control median is 9.597744 min, which maps to nearest zero-based
index 110 on the multiplexed grid. The test allows one output bin in
position and 0.086498 min in retention time, covering the observed
0.086466 min control spread and the compressed scan resolution.
The observed HT apex is index 109 at 9.514768 min, 0.082976 min from
the control median. The archived demultiplexed TIC also peaks at
index 109.
