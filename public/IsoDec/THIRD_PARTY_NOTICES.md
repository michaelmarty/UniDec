# Third-party notices

IsoDec builds and links the IsoGen submodule. IsoGen is distributed under the
BSD 3-Clause License; its license is at `extern/IsoGen/LICENSE`.

IsoGen links FFTW. The bundled FFTW copyright and GPL notices are at
`extern/IsoGen/src/fftw/COPYRIGHT` and `extern/IsoGen/src/fftw/COPYING`.
Distribution wheels bundle FFTW where required by the platform.

Windows release wheels built with Intel oneAPI may include `libmmd.dll` and
`svml_dispmd.dll`. Those files are governed by Intel's runtime redistribution
terms and End User License Agreement.
