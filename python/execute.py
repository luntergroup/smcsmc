""" helper for executing a system command either on the current node or a compute cluster """

import subprocess
import os
import stat

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

    # execute using qrsh
    stdoutname = "{}/{}.stdout".format(dirname, hex(hash(cmd)))
    exitcode = subprocess.check_call("qrsh -j y -cwd -now n -V -N {} -o {} {} {}".format(name, stdoutname, qsub_config, cmd), shell=True)

    # clean up if successful
    if exitcode == 0:
        try:
            os.unlink(stdoutname)
        except OSError:
            pass

    # done
    return exitcode


    
    


