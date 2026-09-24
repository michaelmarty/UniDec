<#
.SYNOPSIS
Synchronizes and pushes the public\IsoDec subtree to michaelmarty/IsoDec main.

.EXAMPLE
.\push_isodec_subtree.ps1

.EXAMPLE
.\push_isodec_subtree.ps1 -WhatIf
#>

[CmdletBinding(SupportsShouldProcess)]
param()

$ErrorActionPreference = "Stop"
$SubtreePrefix = "public/IsoDec"
$Remote = "isodec-public"
$Branch = "main"
$ExpectedRemoteUrl = "https://github.com/michaelmarty/IsoDec.git"

$RepositoryRoot = (& git rev-parse --show-toplevel).Trim()
if ($LASTEXITCODE -ne 0) {
    throw "Run this script from inside the UniDecDev Git repository."
}

Set-Location $RepositoryRoot

if (-not (Test-Path $SubtreePrefix -PathType Container)) {
    throw "IsoDec subtree not found: $SubtreePrefix"
}

$RemoteUrl = (& git remote get-url $Remote).Trim()
if ($LASTEXITCODE -ne 0) {
    throw "Required Git remote '$Remote' is not configured."
}
if ($RemoteUrl -ne $ExpectedRemoteUrl) {
    throw "Remote '$Remote' must point to $ExpectedRemoteUrl; found $RemoteUrl."
}

$UncommittedChanges = & git status --porcelain
if ($LASTEXITCODE -ne 0) {
    throw "Unable to inspect the Git worktree."
}
if ($UncommittedChanges) {
    if ($WhatIfPreference) {
        Write-Warning "A real run requires a clean worktree before synchronizing the subtree."
    }
    else {
        throw "Commit or stash all worktree changes before synchronizing the subtree."
    }
}

if ($PSCmdlet.ShouldProcess("$SubtreePrefix from $Remote/$Branch", "merge remote subtree changes")) {
    & git subtree pull "--prefix=$SubtreePrefix" $Remote $Branch
    if ($LASTEXITCODE -ne 0) {
        throw "Unable to merge $Remote/$Branch into $SubtreePrefix. Resolve any conflicts, commit the merge, and run the script again."
    }
}

if ($PSCmdlet.ShouldProcess("$Remote/$Branch", "push subtree $SubtreePrefix")) {
    & git subtree push --rejoin "--prefix=$SubtreePrefix" $Remote $Branch
    if ($LASTEXITCODE -ne 0) {
        throw "Unable to push $SubtreePrefix to $Remote/$Branch."
    }
}
