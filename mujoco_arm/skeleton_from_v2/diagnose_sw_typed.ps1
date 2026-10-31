$ErrorActionPreference = 'Stop'
$apiPath = 'D:\SolidWorks2026\SOLIDWORKS\api\redist\SolidWorks.Interop.sldworks.dll'
Add-Type -Path $apiPath
Add-Type -ReferencedAssemblies $apiPath -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
using SolidWorks.Interop.sldworks;
public static class SwTypedProbe {
 public static string Probe() {
  object raw=Activator.CreateInstance(Type.GetTypeFromProgID("SldWorks.Application.34",true));
  ISldWorks app=(ISldWorks)raw;
  string version=app.RevisionNumber();
  return "Typed ISldWorks connected; RevisionNumber="+version+"; document count="+app.GetDocumentCount();
 }
}
'@
[SwTypedProbe]::Probe()
