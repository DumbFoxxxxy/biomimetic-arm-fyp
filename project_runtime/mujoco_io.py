"""Load MuJoCo scenes without native Unicode filesystem access.

XML and all file resources are read by Python; native code sees ASCII VFS keys.
Includes retain their document-relative location. No generated mesh is changed.
"""
import hashlib
import xml.etree.ElementTree as ET
from .paths import WorkspacePath as Path, PROJECT_ROOT, project_path


def _candidate(value,base):
    p=Path(value)
    if p.is_absolute():return p
    local=Path(base)/p
    if local.is_file():return local
    # Historical JSON/XML sometimes stored project-relative filenames.
    return project_path(p)


def load_mujoco_xml(xml,assets=None,base_dir=None):
    import mujoco
    if isinstance(xml,bytes):xml=xml.decode('utf-8-sig')
    root=ET.fromstring(xml);base=Path(base_dir or PROJECT_ROOT)
    vfs=dict(assets or {});seen=[]
    compiler=root.find('compiler')
    defaults={} if compiler is None else dict(compiler.attrib)

    def expand(node,document_dir,settings):
        local_compiler=node.find('compiler')
        local_settings=dict(settings)
        if local_compiler is not None:local_settings.update(local_compiler.attrib)
        for child in list(node):
            if child.tag=='include':
                source=_candidate(child.attrib['file'],document_dir)
                if str(source) in seen:raise ValueError('Cyclic MuJoCo include: '+str(source))
                seen.append(str(source));included=ET.fromstring(source.read_text(encoding='utf-8-sig'))
                expand(included,source.parent,local_settings);seen.pop()
                i=list(node).index(child);node.remove(child)
                for item in list(included):node.insert(i,item);i+=1
                continue
            if 'file' in child.attrib:
                raw=child.attrib['file']
                if raw in vfs:continue
                assetdir=local_settings.get('assetdir','')
                kinddir=local_settings.get('meshdir' if child.tag=='mesh' else 'texturedir','')
                folder=kinddir or assetdir
                rawpath=Path(raw)
                source=_candidate(raw,document_dir if rawpath.is_absolute() else Path(document_dir)/folder)
                if not source.is_file():raise FileNotFoundError('MuJoCo resource: '+str(source))
                # Preserve extension for native format detection, use only ASCII names.
                key='asset_'+hashlib.sha256(str(source).encode('utf8')).hexdigest()[:24]+source.suffix.lower()
                if key not in vfs:vfs[key]=source.read_bytes()
                child.set('file',key)
            expand(child,document_dir,local_settings)
        if local_compiler is not None:
            for name in ['assetdir','meshdir','texturedir']:local_compiler.attrib.pop(name,None)

    expand(root,base,defaults)
    return mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'),assets=vfs)


def load_mujoco_model(filename,*args,**kwargs):
    source=project_path(filename)
    return load_mujoco_xml(source.read_text(encoding='utf-8-sig'),base_dir=source.parent,*args,**kwargs)
