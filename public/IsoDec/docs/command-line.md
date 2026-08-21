# Command line

```text
isodec spectrum.txt [--centroided] [--type PEPTIDE|RNA] [-o OUTPUT]
```

When `-o` is omitted, IsoDec writes `<input-stem>_isodec.tsv` beside the input.
Use `--centroided` when the input already contains picked peaks. Complex mass
spectrometry formats require the `unidec` optional dependency.
