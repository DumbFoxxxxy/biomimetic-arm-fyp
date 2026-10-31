"""Scope existing real manufacturing/steel checks to regenerated flexor."""
import check_shaped_muscles as check
check.DEST=check.OUT/'direct_muscle_interfaces/shifted_muscle'
if __name__=='__main__':check.run(muscles=['FJ1486'])
