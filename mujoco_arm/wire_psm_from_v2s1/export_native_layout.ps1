# Strongly typed bridge: creates an independent reference part; no source CAD
# is opened, saved, modified or closed. The point proposal is not hole geometry.
param([switch]$CompileOnly)
$ErrorActionPreference = 'Stop'
$apiPath = 'D:\SolidWorks2026\SOLIDWORKS\api\redist\SolidWorks.Interop.sldworks.dll'
Add-Type -Path $apiPath
Add-Type -ReferencedAssemblies $apiPath -TypeDefinition @'
using System;
using System.IO;
using System.Globalization;
using System.Collections.Generic;
using SolidWorks.Interop.sldworks;
public static class PsmNativeLayout {
 class Marker {public string Module,Role,Owner;public double X,Y,Z;}
 static double Number(string value) {return double.Parse(value,CultureInfo.InvariantCulture)*0.001;}
 static IFeature Finish(IModelDoc2 doc,ISketchManager manager,string name) {
  manager.Insert3DSketch(true);
  IFeature feature=(IFeature)doc.FeatureByPositionReverse(0);
  if(feature==null || !(feature.GetSpecificFeature2() is ISketch))throw new Exception("New sketch missing: "+name);
  feature.Name=name;return feature;
 }
 public static string Run(string folder) {
  var markers=new List<Marker>();var expected=new Dictionary<string,Marker>();
  foreach(string line in File.ReadAllLines(Path.Combine(folder,"anchors_design_mm.csv"))) {
   if(line.StartsWith("module,") || string.IsNullOrWhiteSpace(line))continue;
   string[] c=line.Split(',');markers.Add(new Marker{Module=c[0],Role=c[1],Owner=c[2],X=Number(c[3]),Y=Number(c[4]),Z=Number(c[5])});
  }
  if(markers.Count!=128)throw new Exception("Expected128 marker references");
  ISldWorks app=(ISldWorks)Activator.CreateInstance(Type.GetTypeFromProgID("SldWorks.Application.34",true));
  string template=app.GetDocumentTemplate(1,"",0,0,0);
  if(string.IsNullOrEmpty(template) || !File.Exists(template))throw new Exception("Default SolidWorks part template unavailable");
  IModelDoc2 doc=(IModelDoc2)app.NewDocument(template,0,0,0);
  if(doc==null)throw new Exception("New reference part could not be created");
  ISketchManager manager=doc.SketchManager;
  foreach(Marker p in markers) {
   doc.ClearSelection2(true);manager.Insert3DSketch(true);manager.CreatePoint(p.X,p.Y,p.Z);
   string name="PSM_"+p.Module+"_"+p.Role;Finish(doc,manager,name);expected.Add(name,p);
  }
  string lastGroup=null;double[] previous=null;bool active=false;int lineCount=0;
  foreach(string line in File.ReadAllLines(Path.Combine(folder,"wire_paths_zero_design_mm.csv"))) {
   if(line.StartsWith("module,") || string.IsNullOrWhiteSpace(line))continue;
   string[] c=line.Split(',');string group=c[0]+"_"+c[1];double[] point={Number(c[3]),Number(c[4]),Number(c[5])};
   if(group!=lastGroup) {
    if(active)Finish(doc,manager,"PSM_"+lastGroup+"_ZERO_REFERENCE");
    doc.ClearSelection2(true);manager.Insert3DSketch(true);active=true;previous=null;lastGroup=group;
   }
   if(previous!=null) {
    double dist=Math.Pow(point[0]-previous[0],2)+Math.Pow(point[1]-previous[1],2)+Math.Pow(point[2]-previous[2],2);
    if(dist>1e-20){manager.CreateLine(previous[0],previous[1],previous[2],point[0],point[1],point[2]);lineCount++;}
   }
   previous=point;
  }
  if(active)Finish(doc,manager,"PSM_"+lastGroup+"_ZERO_REFERENCE");
  doc.EditRebuild3();doc.ViewZoomtofit2();string output=Path.Combine(folder,"PSM-W1_REFERENCE_LAYOUT_"+DateTime.Now.ToString("yyyyMMdd_HHmmss")+".SLDPRT");
  int errors=0,warnings=0;bool saved=doc.Extension.SaveAs(output,0,1,null,ref errors,ref warnings);
  if(!saved || errors!=0)throw new Exception("SaveAs failed: "+errors+" warnings="+warnings);
  app.CloseDoc(doc.GetTitle());
  IModelDoc2 check=(IModelDoc2)app.OpenDoc6(output,1,3,"",ref errors,ref warnings);
  if(check==null)throw new Exception("Cannot reopen saved reference part");
  int verified=0;double maximum=0;IFeature f=(IFeature)check.FirstFeature();
  while(f!=null) {
   Marker goal;
   if(expected.TryGetValue(f.Name,out goal)) {
    ISketch sketch=f.GetSpecificFeature2() as ISketch;
    Array pts=sketch==null?null:sketch.GetSketchPoints2() as Array;double nearest=double.PositiveInfinity;
    if(pts!=null)foreach(object element in pts) {
     ISketchPoint p=(ISketchPoint)element;double e=Math.Sqrt(Math.Pow(p.X-goal.X,2)+Math.Pow(p.Y-goal.Y,2)+Math.Pow(p.Z-goal.Z,2));nearest=Math.Min(nearest,e);
    }
    if(nearest>1e-8)throw new Exception("Marker coordinate mismatch: "+f.Name+" "+nearest);
    maximum=Math.Max(maximum,nearest);verified++;
   }
   f=(IFeature)f.GetNextFeature();
  }
  app.CloseDoc(check.GetTitle());
  if(verified!=128)throw new Exception("Expected128 marker sketches, found"+verified);
  string report="{\"native_file\":\""+output.Replace("\\","\\\\")+"\",\"verified_marker_sketches\":"+verified+",\"wire_line_segments\":"+lineCount+",\"maximum_marker_error_m\":"+maximum.ToString("R",CultureInfo.InvariantCulture)+",\"source_CAD_modified\":false,\"reference_layout_only\":true}";
  File.WriteAllText(Path.Combine(folder,"native_export_report.json"),report);return report;
 }
}
'@
if (-not $CompileOnly) { [PsmNativeLayout]::Run((Join-Path (Split-Path -Parent (Split-Path -Parent $PSScriptRoot)) "01_工程骨架方案/输出/wire_psm_from_v2s1/PSM-W1")) }
