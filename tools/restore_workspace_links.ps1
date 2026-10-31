param([switch]$CheckOnly)
$ErrorActionPreference='Stop'
$taskRoot=Split-Path -Parent $PSScriptRoot
$spec=Get-Content -LiteralPath (Join-Path $PSScriptRoot 'workspace_links.json') -Raw | ConvertFrom-Json
foreach ($entry in $spec) {
    $link=[IO.Path]::GetFullPath((Join-Path $taskRoot $entry.path))
    $target=[IO.Path]::GetFullPath((Join-Path $taskRoot $entry.target))
    if (-not $link.StartsWith($taskRoot+'\',[StringComparison]::OrdinalIgnoreCase) -or -not $target.StartsWith($taskRoot+'\',[StringComparison]::OrdinalIgnoreCase)) { throw 'Link outside project root.' }
    if (-not (Test-Path -LiteralPath $target -PathType Container)) {
        Write-Warning "Missing target; skipped: $($entry.target)"
        continue
    }
    if (Test-Path -LiteralPath $link) {
        $item=Get-Item -LiteralPath $link
        if (-not ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) { throw "Existing ordinary directory must be reviewed: $link" }
        $existing=[IO.Path]::GetFullPath([string]$item.Target)
        if ($existing -ne $target) { throw "Existing junction points elsewhere: $link" }
        continue
    }
    if ($CheckOnly) { Write-Output "Would create: $($entry.path)"; continue }
    New-Item -ItemType Directory -Path (Split-Path -Parent $link) -Force | Out-Null
    New-Item -ItemType Junction -Path $link -Target $target | Out-Null
}
Write-Output 'Workspace compatibility link check complete.'
