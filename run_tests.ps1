[CmdletBinding()]
param(
    [string]$Python,
    [string]$ImporterTestData,
    [string]$Hdf5Root = "C:\Program Files\HDF_Group\HDF5\2.1.1",
    [string]$MklRoot = "C:\Program Files (x86)\Intel\oneAPI\mkl\latest",
    [switch]$SkipInstall,
    [switch]$SkipNativeBuild
)

$ErrorActionPreference = "Stop"
$repositoryRoot = $PSScriptRoot
$unidecRoot = Join-Path $repositoryRoot "public\UniDec"

if (-not (Test-Path -LiteralPath $unidecRoot -PathType Container)) {
    throw "UniDec source was not found at '$unidecRoot'."
}

if ([string]::IsNullOrWhiteSpace($Python)) {
    $Python = Join-Path $repositoryRoot ".venv-win\Scripts\python.exe"
}

if ([string]::IsNullOrWhiteSpace($ImporterTestData)) {
    $ImporterTestData = "C:\Python\UniDecImporter\TestData"
}

if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python was not found at '$Python'. Create .venv-win or pass -Python <path>."
}

if (-not (Test-Path -LiteralPath $ImporterTestData -PathType Container)) {
    throw @"
UniDecImporter test data was not found at '$ImporterTestData'.
Clone it into 'C:\Python\UniDecImporter' or pass -ImporterTestData <TestData path>.
"@
}

$fixture = Join-Path $ImporterTestData "IMMS\test_watersimms_txt.txt"
if (-not (Test-Path -LiteralPath $fixture -PathType Leaf) -or (Get-Item -LiteralPath $fixture).Length -eq 0) {
    throw @"
The required IM-MS fixture is missing or empty: '$fixture'.
"@
}

$env:UNIDEC_IMPORTER_TEST_DATA = (Resolve-Path -LiteralPath $ImporterTestData).Path
$env:UNIDEC_REQUIRE_IMPORTER_TEST_DATA = "1"

Push-Location $unidecRoot
try {
    if (-not $SkipNativeBuild) {
        $hdf5RootPath = (Resolve-Path -LiteralPath $Hdf5Root -ErrorAction Stop).Path
        $mklRootPath = (Resolve-Path -LiteralPath $MklRoot -ErrorAction Stop).Path
        $intelSetup = "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat"
        $intelCompilerSetup = "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat"
        $mklConfig = Join-Path $mklRootPath "lib\cmake\mkl\MKLConfig.cmake"
        foreach ($file in @($intelSetup, $intelCompilerSetup, $mklConfig)) {
            if (-not (Test-Path -LiteralPath $file -PathType Leaf)) {
                throw "Native build prerequisite was not found: '$file'."
            }
        }
        foreach ($command in @("cmake.exe", "ninja.exe", "wsl.exe")) {
            if (-not (Get-Command $command -ErrorAction SilentlyContinue)) {
                throw "Native build prerequisite was not found on PATH: $command."
            }
        }

        $runId = [guid]::NewGuid().ToString("N").Substring(0, 8)
        $windowsBuild = Join-Path $unidecRoot "build\local-intel-$runId"
        New-Item -ItemType Directory -Path $windowsBuild -Force | Out-Null
        $hdf5Path = $hdf5RootPath.Replace('\', '/')
        $mklPath = $mklRootPath.Replace('\', '/')
        $sourcePath = (Join-Path $unidecRoot "unidec\src").Replace('\', '/')
        $buildPath = $windowsBuild.Replace('\', '/')

        $buildCommand = Join-Path $windowsBuild "build-intel.cmd"
        @"
@echo off
call "$intelSetup" vs2022
if errorlevel 1 exit /b 1
call "$intelCompilerSetup"
if errorlevel 1 exit /b 1
set "CC=icx-cl"
set "CXX=icx-cl"
cmake -S "$sourcePath" -B "$buildPath" -G Ninja -DCMAKE_BUILD_TYPE=Release -DUNIDEC_REQUIRE_INTELLLVM=ON -DHDF5_ROOT="$hdf5Path" -DMKL_ROOT="$mklPath" -DMKL_DIR="$mklPath/lib/cmake/mkl"
if errorlevel 1 exit /b 1
cmake --build "$buildPath" --parallel
if errorlevel 1 exit /b 1
"@ | Set-Content -LiteralPath $buildCommand -Encoding Ascii
        & $env:ComSpec /d /c $buildCommand
        if ($LASTEXITCODE -ne 0) { throw "The IntelLLVM Windows build failed." }

        $nativeBin = Join-Path $unidecRoot "unidec\bin"
        foreach ($name in @("unidec.exe", "unideclib.dll")) {
            $built = Join-Path $windowsBuild $name
            $deployed = Join-Path $nativeBin $name
            if (-not (Test-Path -LiteralPath $built -PathType Leaf) -or
                -not (Test-Path -LiteralPath $deployed -PathType Leaf) -or
                (Get-FileHash -LiteralPath $built).Hash -ne (Get-FileHash -LiteralPath $deployed).Hash) {
                throw "The current Windows build was not deployed to '$deployed'."
            }
        }

        $runtimeDlls = foreach ($source in @(
            (Join-Path $hdf5RootPath "bin\hdf5.dll")
            (Join-Path $hdf5RootPath "bin\hdf5_hl.dll")
            (Join-Path $hdf5RootPath "bin\vcruntime140.dll")
            (Join-Path $mklRootPath "bin\mkl_rt.2.dll")
            (Join-Path $mklRootPath "bin\mkl_core.2.dll")
            (Join-Path $mklRootPath "bin\mkl_intel_thread.2.dll")
            (Join-Path $mklRootPath "bin\mkl_def.2.dll")
            (Join-Path $mklRootPath "bin\mkl_avx2.2.dll")
            "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libiomp5md.dll"
            "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\libmmd.dll"
            "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\svml_dispmd.dll"
        )) {
            Get-Item -LiteralPath $source -ErrorAction Stop
        }
        Copy-Item -LiteralPath $runtimeDlls.FullName -Destination $nativeBin -Force

        $savedPath = $env:PATH
        try {
            $env:PATH = "$nativeBin;$env:SystemRoot\System32;$env:SystemRoot"
            $nativeHelp = & (Join-Path $nativeBin "unidec.exe") 2>&1
            if ($LASTEXITCODE -ne 88 -or ($nativeHelp -join "`n") -notmatch "Universal Deconvolution") {
                throw "The Windows native executable could not load using only its bundled runtime DLLs."
            }
        }
        finally {
            $env:PATH = $savedPath
        }

        $wslRoot = & wsl.exe --exec wslpath -a $unidecRoot
        if ($LASTEXITCODE -ne 0 -or -not $wslRoot) { throw "WSL could not locate the UniDec checkout." }
        $wslRoot = $wslRoot.Trim()
        $linuxBuild = "$wslRoot/build/local-wsl-$runId"
        & wsl.exe --exec env "UNIDEC_BUILD_DIR=$linuxBuild" bash "$wslRoot/build_native.sh"
        if ($LASTEXITCODE -ne 0) { throw "The WSL Linux native build failed." }
        foreach ($name in @("unideclinux", "libunideclib.so")) {
            & wsl.exe --exec cmp "$linuxBuild/$name" "$wslRoot/unidec/bin/$name"
            if ($LASTEXITCODE -ne 0) { throw "The current WSL build was not deployed: $name." }
            $linuxDependencies = & wsl.exe --exec ldd "$linuxBuild/$name" 2>&1
            if ($LASTEXITCODE -ne 0 -or ($linuxDependencies -join "`n") -match "not found") {
                throw "The WSL native output has an unresolved shared library: $name.`n$($linuxDependencies -join "`n")"
            }
        }
        Write-Host "Current Windows and Linux native binaries are bundled in unidec/bin."
    }

    if (-not $SkipInstall) {
        & $Python -m pip install --upgrade pip
        if ($LASTEXITCODE -ne 0) {
            exit $LASTEXITCODE
        }

        & $Python -m pip install -e ".[dev]"
        if ($LASTEXITCODE -ne 0) {
            exit $LASTEXITCODE
        }
    }

    & $Python -m unittest discover -s tests -p "test_*.py" -v
    exit $LASTEXITCODE
}
finally {
    Pop-Location
}
