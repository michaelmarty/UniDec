:orphan:

UniDec archive index
====================

``native-src-2026-11-25.zip``
  Snapshot of the native archive directory. Contains four historical C
  headers (``UD_match.h``, ``UD_peak_width.h``, ``UD_sg.h``, and
  ``UniDecLC_Main.h``), a UniChrom performance plan, and notes and a patch
  for restoring the retired UCCD-style blur comparison mode.

  The ZIP retains the original UniDecDev paths under
  ``public/UniDec/unidec/src/Archive/``. Extract it outside your working
  tree and read ``UNICHROM_UCCD_BLUR.md`` before trying the patch.

``legacy-visual-studio-project-2026-11-25.zip``
  The former ``UniDec.sln`` and ``UniDec.vcxproj``, together at their
  original paths under ``public/UniDec/unidec/src/``.
  Use ``unidec/src/CMakeLists.txt`` and the native build instructions in
  the UniDec README for current builds.
