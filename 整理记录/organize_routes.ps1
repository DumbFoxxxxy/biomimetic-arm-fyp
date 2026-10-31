$ErrorActionPreference = 'Stop'
$taskRoot = 'D:\XJTLU\PSP'
$oldRoute = Join-Path $taskRoot '01_工程骨架方案'
$newRoute = Join-Path $taskRoot '02_人体骨骼晶格肌肉方案'
$common = Join-Path $taskRoot '00_共用资料'
$recordPath = Join-Path $taskRoot '整理记录\migration.json'
if (Test-Path -LiteralPath $recordPath) {
    $prior = @(Get-Content -LiteralPath $recordPath -Raw | ConvertFrom-Json)
    foreach ($r in $prior) {
        if ($r.state -ne 'prepared' -or -not (Test-Path -LiteralPath $r.source) -or (Test-Path -LiteralPath $r.destination)) { throw 'Existing migration needs manual state review.' }
        if ((Get-Item -LiteralPath $r.source).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Unexpected source link during resume.' }
    }
    Copy-Item -LiteralPath $recordPath -Destination (Join-Path $taskRoot '整理记录\migration_rolled_back_attempt.json')
}
foreach ($p in @($oldRoute,$newRoute,$common)) { New-Item -ItemType Directory -Path $p -Force | Out-Null }
$moves = @(
    @('final_model','01_工程骨架方案\SolidWorks'),
    @('model','01_工程骨架方案\早期模型'),
    @('archive','01_工程骨架方案\历史归档'),
    @('骨骼','02_人体骨骼晶格肌肉方案\输入\骨骼'),
    @('muscle','02_人体骨骼晶格肌肉方案\输入\muscle')
)
foreach ($dir in Get-ChildItem -LiteralPath (Join-Path $taskRoot 'mujoco_arm\exports') -Directory) {
    $branch = if ($dir.Name -in @('anatomical_arm','anatomical_concept')) { '02_人体骨骼晶格肌肉方案' } else { '01_工程骨架方案' }
    $moves += ,@(('mujoco_arm\exports\'+$dir.Name),($branch+'\输出\'+$dir.Name))
}
$records = @()
foreach ($move in $moves) {
    $src = Join-Path $taskRoot $move[0]; $dst = Join-Path $taskRoot $move[1]
    $resolved = (Resolve-Path -LiteralPath $src).Path
    if (-not $resolved.StartsWith($taskRoot+'\',[StringComparison]::OrdinalIgnoreCase)) { throw "Source outside workspace: $resolved" }
    $fullDestination = [IO.Path]::GetFullPath($dst)
    if (-not $fullDestination.StartsWith($taskRoot+'\',[StringComparison]::OrdinalIgnoreCase)) { throw "Destination outside workspace: $fullDestination" }
    if ((Get-Item -LiteralPath $src).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw "Unexpected existing link: $src" }
    if (Test-Path -LiteralPath $dst) { throw "Destination exists: $dst" }
    # Get-ChildItem does not follow junctions unless explicitly requested.
    $files = @(Get-ChildItem -LiteralPath $src -File -Recurse)
    $inventory = @($files | ForEach-Object { @{relative=$_.FullName.Substring($src.Length+1); length=$_.Length; modified_ticks=$_.LastWriteTimeUtc.Ticks} })
    $record = @{source=$src; destination=$dst; files=$inventory; count=$files.Count; state='prepared'}
    $records += $record
    $records | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $recordPath -Encoding utf8
    New-Item -ItemType Directory -Path (Split-Path -Parent $dst) -Force | Out-Null
    Move-Item -LiteralPath $src -Destination $dst
    try { New-Item -ItemType Junction -Path $src -Target $dst | Out-Null }
    catch { Move-Item -LiteralPath $dst -Destination $src; throw }
    foreach ($f in $inventory) {
        $item = Get-Item -LiteralPath (Join-Path $dst $f.relative)
        if ($item.Length -ne $f.length -or $item.LastWriteTimeUtc.Ticks -ne $f.modified_ticks) { throw "Inventory mismatch: $($item.FullName)" }
    }
    $record.state='moved_verified_compatibility_junction_present'
    $records | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $recordPath -Encoding utf8
    Write-Output ("Moved {0}: {1} files" -f $move[0],$files.Count)
}
foreach ($dir in Get-ChildItem -LiteralPath (Join-Path $taskRoot 'mujoco_arm') -Directory) {
    if ($dir.Name -in @('exports','__pycache__')) { continue }
    $branch = if ($dir.Name -in @('anatomical_arm','anatomical_concept')) { $newRoute } else { $oldRoute }
    $link = Join-Path $branch ('脚本\'+$dir.Name)
    New-Item -ItemType Directory -Path (Split-Path -Parent $link) -Force | Out-Null
    New-Item -ItemType Junction -Path $link -Target $dir.FullName | Out-Null
}
foreach ($name in @('knowledge_base','literature_inbox','design','tmp')) {
    New-Item -ItemType Junction -Path (Join-Path $common $name) -Target (Join-Path $taskRoot $name) | Out-Null
}
New-Item -ItemType Junction -Path (Join-Path $oldRoute '旧演示代码入口') -Target (Join-Path $taskRoot 'mujoco_arm') | Out-Null
Write-Output 'Route data relocated; original paths retained as compatibility junctions. Source code unchanged.'
