"""Standalone bench interfaces of the revised flexor, not installed-arm proof."""
from build_direct_muscle_yoke import *
import check_muscle_end_assembly as ends
if __name__=='__main__':
    fid='FJ1477' if EXTENSOR_REANCHOR else 'FJ1486'
    ends.DEST=muscle_folder(fid).parent
    ends.run((fid,))
