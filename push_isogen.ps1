[CmdletBinding()]
param(
    [Parameter(Mandatory)]
    [ValidateNotNullOrEmpty()]
    [string]$Message,

    [string]$Branch
)

$ErrorActionPreference = "Stop"

function Invoke-Git {
    param([string[]]$Arguments)

    Write-Host "+ git -C $IsoGenPath $($Arguments -join ' ')"
    & git -C $IsoGenPath @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Git command failed: git -C $IsoGenPath $($Arguments -join ' ')"
    }
}

$IsoGenPath = Join-Path $PSScriptRoot "extern\IsoGen"
if (-not (Test-Path (Join-Path $IsoGenPath ".git"))) {
    throw "IsoGen submodule is not initialized. Run: git submodule update --init --recursive"
}

Invoke-Git @("fetch", "origin", "--prune")

if (-not $Branch) {
    $defaultRef = (& git -C $IsoGenPath symbolic-ref --quiet --short refs/remotes/origin/HEAD).Trim()
    if ($LASTEXITCODE -ne 0 -or -not $defaultRef.StartsWith("origin/")) {
        throw "Could not determine IsoGen's default branch. Specify it with -Branch."
    }
    $Branch = $defaultRef.Substring("origin/".Length)
}

& git -C $IsoGenPath show-ref --verify --quiet "refs/remotes/origin/$Branch"
if ($LASTEXITCODE -ne 0) {
    throw "Remote branch origin/$Branch does not exist."
}

& git -C $IsoGenPath show-ref --verify --quiet "refs/heads/$Branch"
if ($LASTEXITCODE -eq 0) {
    Invoke-Git @("switch", $Branch)
}
else {
    Invoke-Git @("switch", "--track", "-c", $Branch, "origin/$Branch")
}

Invoke-Git @("pull", "--ff-only", "origin", $Branch)

$changes = (& git -C $IsoGenPath status --porcelain | Out-String).Trim()
if (-not $changes) {
    Write-Host "IsoGen has no local changes to commit; nothing was pushed."
    Write-Host "If this release updates IsoDec's submodule pointer, commit extern/IsoGen from the repository root."
    return
}

Invoke-Git @("add", "--all")
Invoke-Git @("commit", "-m", $Message)
Invoke-Git @("push", "origin", $Branch)

Write-Host "Pushed IsoGen commit $((& git -C $IsoGenPath rev-parse --short HEAD).Trim()) to origin/$Branch."
