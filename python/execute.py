""" helper for executing a system command either on the current node or a compute cluster """

import subprocess
import os

# global configuration settings; None to not submit jobs to cluster
qsub_config = None


def check_call(cmd, shell=True, outputdir=None, name=None):
    if qsub_config is None:
        return subprocess.check_call(cmd, shell=True)

    if outputdir is None: outputdir = "./"
    if name is None: name = cmd.split(' ')[0]

    # create output directory if necessary
    dirname = os.path.abspath( outputdir )
    if not os.path.isdir( dirname ):
        try:
            os.makedirs( dirname )
        except:
            raise ValueError("Cannot find or create directory '{}'".format( dirname ))

    # create shell script in output directory
    scriptname = "{}/run.{}.sh".format(dirname, hex(hash(cmd)))
    stdoutname = "{}/{}.stdout".format(dirname, hex(hash(cmd)))
    script = open(scriptname,'w')
    script.write("#!/bin/bash\n")
    script.write("#$ -j y\n")
    script.write("#$ -o {}\n".format(stdoutname))
    script.write("#$ -sync y\n")
    script.write("#$ -cwd\n")
    script.write("#$ -V\n")
    script.write("#$ -N {}\n".format(name))
    script.write("#$ {}\n\n".format(qsub_config))
    script.write("{}\n".format(cmd))
    script.close()
    
    # submit script, wait until completion, and obtain exit code (because we set -sync y)
    exitcode = subprocess.check_call("qsub {}".format(scriptname), shell=True)

    # clean up if successful
    if exitcode == 0:
        os.unlink(scriptname)
        try:
            os.unlink(stdoutname)
        except OSError:
            pass

    # done
    return exitcode


    
    


