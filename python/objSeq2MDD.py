import meshio
import glob
import numpy as np
from struct import pack

dir = "/mnt/F/jjcasmar/projects/VNCS_Scenes/AnimatedCylinder_SCA_6/Coarse/Simulation/"

fines = glob.glob(dir + "fine_**.vtu.vtu")
if len(fines):
    fine0 = meshio.read(fines[0])
    fine0.write(dir + "fine.obj")

    numverts = len(fine0.points)
    numframes = len(fines)
    print(numframes)

    f = open(dir + "/seq_offsets.mdd", 'wb')
    f.write(pack(">2i", numframes, numverts))
    f.write(pack(">%df" % numframes, *[frame / (240.0)
            for frame in range(numframes)]))

    f.write(pack(">%df" % (numverts * 3), *
            [axis for v in fine0.points for axis in v]))

    for fine in sorted(fines):
        fineMesh = meshio.read(fine)
        f.write(pack(">%df" % (numverts * 3), *
                [axis for v in fineMesh.points for axis in v]))
