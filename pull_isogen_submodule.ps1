<#
.SYNOPSIS
Updates the IsoGen submodule and attaches it to its main branch.

.EXAMPLE
.\pull_isogen_submodule.ps1
#>

[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$SubmodulePath = "public/IsoDec/extern/IsoGen"
$Branch = "main"
$IsoGenPath = Join-Path $PSScriptRoot $SubmodulePath

function Invoke-IsoGenGit {
    param([string[]]$Arguments)

    Write-Host "+ git -C $IsoGenPath $($Arguments -join ' ')"
    & git -C $IsoGenPath @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Git command failed: git -C $IsoGenPath $($Arguments -join ' ')"
    }
}

if (-not (Test-Path -LiteralPath (Join-Path $PSScriptRoot ".git"))) {
    throw "The script must remain in the UniDecDev repository root."
}

if (-not (Test-Path -LiteralPath (Join-Path $IsoGenPath ".git"))) {
    Write-Host "+ git submodule update --init --recursive -- $SubmodulePath"
    & git -C $PSScriptRoot submodule update --init --recursive -- $SubmodulePath
    if ($LASTEXITCODE -ne 0) {
        throw "Unable to initialize the IsoGen submodule."
    }
}

$Changes = (& git -C $IsoGenPath status --porcelain | Out-String).Trim()
if ($LASTEXITCODE -ne 0) {
    throw "Unable to inspect the IsoGen worktree."
}
if ($Changes) {
    throw "IsoGen has uncommitted changes. Commit or stash them before updating."
}

Invoke-IsoGenGit @("fetch", "origin", "--prune", "--tags")

& git -C $IsoGenPath show-ref --verify --quiet "refs/remotes/origin/$Branch"
if ($LASTEXITCODE -ne 0) {
    throw "Remote branch origin/$Branch does not exist."
}

& git -C $IsoGenPath symbolic-ref --quiet HEAD *> $null
if ($LASTEXITCODE -ne 0) {
    & git -C $IsoGenPath merge-base --is-ancestor HEAD "origin/$Branch"
    if ($LASTEXITCODE -ne 0) {
        throw "Detached HEAD contains a commit outside origin/$Branch. Create a branch for it before updating."
    }
}

& git -C $IsoGenPath show-ref --verify --quiet "refs/heads/$Branch"
if ($LASTEXITCODE -eq 0) {
    Invoke-IsoGenGit @("switch", $Branch)
    Invoke-IsoGenGit @("branch", "--set-upstream-to=origin/$Branch", $Branch)
}
else {
    Invoke-IsoGenGit @("switch", "--track", "-c", $Branch, "origin/$Branch")
}

Invoke-IsoGenGit @("pull", "--ff-only", "origin", $Branch)

$Version = (& git -C $IsoGenPath describe --tags --exact-match HEAD 2>$null | Out-String).Trim()
if (-not $Version) {
    $Version = (& git -C $IsoGenPath rev-parse --short HEAD).Trim()
}

Write-Host "IsoGen is on $Branch at $Version and tracks origin/$Branch."
Write-Host "If the submodule pointer changed, commit $SubmodulePath in UniDecDev."
