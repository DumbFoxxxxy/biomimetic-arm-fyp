# Uses IDispatch directly: independent of a broken installed type library.
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
using DISPPARAMS=System.Runtime.InteropServices.ComTypes.DISPPARAMS;
using EXCEPINFO=System.Runtime.InteropServices.ComTypes.EXCEPINFO;
public static class CadDispatch {
 [ComImport, Guid("00020400-0000-0000-C000-000000000046"), InterfaceType(ComInterfaceType.InterfaceIsIUnknown)]
 interface Dispatch {
  [PreserveSig] int GetTypeInfoCount(out uint count);
  [PreserveSig] int GetTypeInfo(uint index,uint lcid,out IntPtr info);
  [PreserveSig] int GetIDsOfNames(ref Guid iid,[MarshalAs(UnmanagedType.LPArray,ArraySubType=UnmanagedType.LPWStr)] string[] names,uint count,uint lcid,[Out] int[] ids);
  [PreserveSig] int Invoke(int id,ref Guid iid,uint lcid,ushort flags,ref DISPPARAMS args,[MarshalAs(UnmanagedType.Struct)] out object result,ref EXCEPINFO exception,out uint argError);
 }
 [DllImport("oleaut32.dll")] static extern int VariantClear(IntPtr ptr);
 public static object Call(object instance,string name,object[] args) {
  var disp=(Dispatch)instance;var iid=Guid.Empty;var ids=new int[1];
  Marshal.ThrowExceptionForHR(disp.GetIDsOfNames(ref iid,new[]{name},1,1033,ids));
  int stride=IntPtr.Size==8?24:16;var mem=Marshal.AllocCoTaskMem(stride*args.Length);
  var p=new DISPPARAMS{rgvarg=mem,cArgs=args.Length};
  try {
   for(int i=0;i<args.Length;i++)Marshal.GetNativeVariantForObject(args[args.Length-i-1],IntPtr.Add(mem,i*stride));
   object result;var ex=new EXCEPINFO();uint bad;
   int hr=disp.Invoke(ids[0],ref iid,1033,1,ref p,out result,ref ex,out bad);
   if(hr<0)throw new Exception(name+": "+ex.bstrDescription+" HRESULT="+hr.ToString("X"));
   return result;
  } finally {for(int i=0;i<args.Length;i++)VariantClear(IntPtr.Add(mem,i*stride));Marshal.FreeCoTaskMem(mem);}
 }
}
'@
$ErrorActionPreference = 'Stop'
foreach ($cadProgId in @('SldWorks.Application.34','SldWorks.Application')) {
 try {
  $swRaw = [Runtime.InteropServices.Marshal]::GetActiveObject($cadProgId)
  Write-Output "$cadProgId connected"
  [CadDispatch]::Call($swRaw,'RevisionNumber',@())
 } catch { Write-Output "$cadProgId : $($_.Exception.Message)" }
}
Get-ItemProperty 'Registry::HKEY_CLASSES_ROOT/SldWorks.Application/CLSID','Registry::HKEY_CLASSES_ROOT/SldWorks.Application/CurVer','Registry::HKEY_CLASSES_ROOT/SldWorks.Application.34/CLSID' -ErrorAction SilentlyContinue | Select-Object '(default)',PSPath
