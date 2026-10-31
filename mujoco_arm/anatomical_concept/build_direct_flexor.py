"""Generate the actual shifted TPU lattice/shell/guides/pads and cap hardware."""
import build_shaped_elbow_muscles as generator
generator.DEST=generator.OUT/'direct_muscle_interfaces/shifted_muscle'
if __name__=='__main__':generator.run(muscles=['FJ1486'])
