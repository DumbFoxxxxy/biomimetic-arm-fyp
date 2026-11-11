"""Run the common spool driver with explicitly recorded moving-bone contact forces."""
import argparse
import sys
from build_s1_teres_origin_mount import ROOT,read,save
from build_channel_first_inventory import provenance
from channel_first_moving_sheath_contact import MovingSheathContact
import simulate_channel_first_elbow_channel_beams as driver


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',required=True)
    parser.add_argument('--wall-contact-density',type=float,default=200000.)
    parser.add_argument('--wall-friction',type=float,default=0.)
    parser.add_argument('--wall-update-interval',type=float,default=.0002)
    args,remaining=parser.parse_known_args()
    assert args.wall_contact_density>0 and args.wall_friction>=0 and args.wall_update_interval>0
    check_file=ROOT/args.folder/'model_check.json';check=read(check_file)
    assert check.get('external_clamps') and check.get('wire_paths')
    driver_source=provenance(ROOT/'tools/simulate_channel_first_elbow_channel_beams.py')
    contact=[None];original_pad=driver.pad_forces;original_save=driver.save
    def augmented_support(model,data,pads,scale,damping):
        loads=original_pad(model,data,pads,scale,damping)
        if contact[0] is None:
            contact[0]=MovingSheathContact(model,check,args.wall_contact_density,args.wall_friction,args.wall_update_interval)
        contact[0].apply(data)
        return loads
    def save_with_contact(path,value):
        # Record the wrapper and integrated driver, so extra forces cannot be
        #mistaken for the no-contact protocol when replaying the same XML.
        value['integration_driver']=driver_source
        value['generator']=provenance(ROOT/'tools/simulate_channel_first_external_contact.py')
        value['moving_sheath_bone_contact']=contact[0].report() if contact[0] else None
        value['contact_protocol_note']='Distributed external sheath forces added after scapular pads at each step. Held generalized forces between recorded updates. XML alone does not apply these forces. Bone friction is not internal cable friction.'
        original_save(path,value)
    driver.pad_forces=augmented_support;driver.save=save_with_contact
    sys.argv=['simulate_channel_first_elbow_channel_beams','--folder',args.folder,*remaining]
    try:driver.main()
    except Exception as error:
        save(check_file.parent/'external_contact_failure.json',dict(model_check=provenance(check_file),
            generator=provenance(ROOT/'tools/simulate_channel_first_external_contact.py'),integration_driver=driver_source,
            error=repr(error),last_contact_time_s=contact[0].last_time if contact[0] else None,
            contact_diagnostics=contact[0].report() if contact[0] else None,
            status='failed_run_not_motion_acceptance'))
        raise


if __name__=='__main__':main()
