$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$partDir = Join-Path $root "01_工程骨架方案\SolidWorks"
$rawDir = Join-Path $root "mujoco_arm\assets_tendon_demo\raw"
$boxPath = Join-Path $root "mujoco_arm\assets_tendon_demo\body_boxes.json"
New-Item -ItemType Directory -Force -Path $rawDir | Out-Null

$code = @'
using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

public static class SwTendonDemoExporter {
    static object Call(object o, string name, params object[] args) {
        return o.GetType().InvokeMember(
            name, BindingFlags.InvokeMethod, null, o, args
        );
    }

    static string ArrayText(Array values) {
        var parts = new List<string>();
        foreach (object value in values) {
            parts.Add(Convert.ToString(value));
        }
        return string.Join(",", parts);
    }

    public static string Export(string partPath, string stlPath) {
        Type type = Type.GetTypeFromProgID("SldWorks.Application");
        object sw = Activator.CreateInstance(type);
        object doc = Call(sw, "OpenDoc", partPath, 1);
        if (doc == null) {
            throw new Exception("Cannot open " + partPath);
        }
        try {
            Array bodies = Call(doc, "GetBodies2", 0, true) as Array;
            if (bodies == null || bodies.Length == 0) {
                throw new Exception("No solid body in " + partPath);
            }
            Array box = Call(bodies.GetValue(0), "GetBodyBox") as Array;
            object saved = Call(doc, "SaveAs", stlPath);
            if (!Convert.ToBoolean(saved)) {
                throw new Exception("STL export failed for " + partPath);
            }
            return ArrayText(box);
        }
        finally {
            Call(sw, "CloseDoc", Call(doc, "GetTitle"));
        }
    }
}
'@

Add-Type -TypeDefinition $code -Language CSharp

$rows = @()
foreach ($name in @("thorax", "clavicle", "scapula", "humerus", "ulna")) {
    $partPath = Join-Path $partDir "$name.SLDPRT"
    if (-not (Test-Path -LiteralPath $partPath)) {
        $partPath = Join-Path $partDir ($name.Substring(0, 1).ToUpper() + $name.Substring(1) + ".SLDPRT")
    }
    if (-not (Test-Path -LiteralPath $partPath)) {
        throw "Missing part: $name"
    }
    $stlPath = Join-Path $rawDir "$name.stl"
    $boxText = [SwTendonDemoExporter]::Export($partPath, $stlPath)
    $values = $boxText.Split(",") | ForEach-Object { [double]$_ }
    $rows += [pscustomobject]@{
        name = $name
        box_min_m = @($values[0], $values[1], $values[2])
        box_max_m = @($values[3], $values[4], $values[5])
        stl = $stlPath
    }
    Write-Output "exported $name"
}

$rows | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $boxPath -Encoding UTF8
Write-Output "wrote $boxPath"
