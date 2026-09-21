[CmdletBinding()]
param(
    [string]$Python,
    [string]$ImporterTestData,
    [switch]$SkipInstall
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
