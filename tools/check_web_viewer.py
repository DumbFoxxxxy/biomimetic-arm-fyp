"""Exercise the real local WebGL page and record desktop/mobile screenshots."""
from functools import partial
import hashlib
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
import json
from pathlib import Path
import threading
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'reviews/S1_2026-11-03/web_viewer'

class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self,*args):pass

def main():
    DEST.mkdir(parents=True,exist_ok=True)
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(QuietHandler,directory=str(ROOT/'docs')))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    errors=[];bad_responses=[];checks={}
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(executable_path='C:/Program Files/Google/Chrome/Application/chrome.exe',headless=True,
                args=['--use-angle=swiftshader','--enable-unsafe-swiftshader'])
            page=browser.new_page(viewport={'width':1440,'height':960},device_scale_factor=1)
            page.on('pageerror',lambda error:errors.append(str(error)))
            page.on('response',lambda response:bad_responses.append(dict(url=response.url,status=response.status)) if response.status>=400 else None)
            page.goto(f'http://127.0.0.1:{server.server_port}/',wait_until='networkidle',timeout=90000)
            page.wait_for_function('window.__armViewer?.ready === true',timeout=90000)
            page.wait_for_function('window.__armViewer.renderer.info.render.frame > 2')
            get=lambda expression:page.evaluate('() => {const v=window.__armViewer;return '+expression+';}')
            state=get('({parts:v.manifest.part_count,triangles:v.manifest.triangle_count,rendered:v.renderer.info.render.triangles,meshes:Object.values(v.layers).reduce((s,g)=>s+g.children.length,0)})')
            assert state['parts']==state['meshes']==96 and state['rendered']==state['triangles']==682744,state
            keys=set(get('v.manifest.parts.map(p=>p.key)'))
            assert {'bone:humerus','bone:forearm','tpu:elbow_flexor','tpu:elbow_extensor'} <= keys
            checks['full_scene']=state
            before=get('v.camera.position.toArray()');box=page.locator('canvas').bounding_box()
            x=box['x']+box['width']*.5;y=box['y']+box['height']*.45
            page.mouse.move(x,y);page.mouse.down();page.mouse.move(x+130,y+45,steps=12);page.mouse.up();page.wait_for_timeout(400)
            after=get('v.camera.position.toArray()');assert max(abs(a-b) for a,b in zip(before,after))>1e-4
            checks['mouse_orbit_changes_camera']=True
            distance=get('v.camera.position.distanceTo(v.controls.target)');page.mouse.wheel(0,-300);page.wait_for_timeout(400)
            assert get('v.camera.position.distanceTo(v.controls.target)')<distance
            checks['wheel_zooms']=True
            before=get('v.controls.target.toArray()');page.mouse.down(button='right');page.mouse.move(x+190,y+90,steps=10);page.mouse.up(button='right');page.wait_for_timeout(400)
            assert max(abs(a-b) for a,b in zip(before,get('v.controls.target.toArray()')))>1e-5
            checks['right_drag_pans']=True
            page.get_by_role('button',name='骨骼结构',exact=True).click()
            visibility=get('Object.fromEntries(Object.entries(v.layers).map(([k,g])=>[k,g.visible]))')
            assert visibility==dict(bones=True,muscles=False,wires=False,hardware=False)
            checks['skeleton_preset']=visibility
            page.locator('[data-layer="muscles"]').check();assert get('v.layers.muscles.visible')
            page.get_by_role('button',name='肌肉布局',exact=True).click();assert not get('v.layers.hardware.visible') and get('v.layers.wires.visible')
            checks['layer_toggles_and_muscle_preset']=True
            page.get_by_role('button',name='自动旋转',exact=True).click();assert get('v.controls.autoRotate')
            before=get('v.camera.position.toArray()');page.wait_for_timeout(700)
            assert max(abs(a-b) for a,b in zip(before,get('v.camera.position.toArray()')))>1e-5
            page.get_by_role('button',name='自动旋转',exact=True).click();checks['auto_rotate']=True
            page.get_by_role('button',name='完整装配',exact=True).click();page.get_by_role('button',name='重置视角',exact=True).click();page.wait_for_timeout(400)
            assert max(abs(x) for x in get('v.controls.target.toArray()'))<1e-8
            checks['reset_view']=True
            page.screenshot(path=str(ROOT/'docs/viewer/preview.png'),full_page=True)
            if not page.get_by_role('button',name='全屏',exact=True).is_disabled():
                page.get_by_role('button',name='全屏',exact=True).click();page.wait_for_timeout(200)
                assert page.evaluate('document.fullscreenElement?.id')=='viewport'
                page.get_by_role('button',name='全屏',exact=True).click();page.wait_for_timeout(200);checks['fullscreen']=True
            page.set_viewport_size({'width':390,'height':844});page.wait_for_timeout(500);page.get_by_role('button',name='重置视角',exact=True).click();page.wait_for_timeout(400)
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
            assert page.locator('canvas').bounding_box()['width']>=300
            checks['mobile_no_horizontal_overflow']=True
            page.screenshot(path=str(DEST/'mobile.png'),full_page=True)
            page.locator('#shoulder-view').click()
            page.wait_for_function('window.__armViewer?.ready === true',timeout=90000)
            assert get('v.manifest.part_count')==191
            assert 'bone:scapula' in get('v.manifest.parts.map(p=>p.key)')
            assert 'bone:forearm' not in get('v.manifest.parts.map(p=>p.key)')
            checks['separate_shoulder_scene']=True
            assert not errors and not bad_responses,(errors,bad_responses)
            browser.close()
    finally:server.shutdown();server.server_close()
    result=dict(scope='Real Chrome WebGL page and interaction checks; no physical actuation or motion simulation.',checks=checks,
        page_errors=errors,bad_asset_responses=bad_responses,
        desktop_preview='docs/viewer/preview.png',mobile_preview='reviews/S1_2026-11-03/web_viewer/mobile.png')
    (DEST/'browser_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,ensure_ascii=True,indent=2),flush=True)

if __name__=='__main__':main()
