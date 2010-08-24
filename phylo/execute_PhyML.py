# this is the basic template for executing a particular command on all the files in a directory.
# I just modify as necessary and run it at the python prompt (=ipython)

import glob, shlex, subprocess
paths = glob.glob('path/to/files/*.phylip')
for p in paths:
    command = "bsub -q normal_serial -R 'select[mem>3000]' 'PhyML_3.0_linux64 -i %s -m GTR'" % (p)
    command = shlex.split(command)
    subprocess.Popen(command)