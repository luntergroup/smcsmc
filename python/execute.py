""" helper for executing a system command either on the current node or a compute cluster """

import subprocess
import os
import stat

# global configuration settings; None to not submit jobs to cluster
qsub_config = None


def check_call(cmd, shell=True, outputdir=None, name=None, use_popen=False, use_submit=False, submit_local=False):

    if submit_local or (qsub_config is None):
        return subprocess.check_call(cmd, shell=True)

    if outputdir is None: outputdir = "./"
    if name is None: name = cmd.split(' ')[0]
    name = name.split('/')[-1]

    # create output directory if necessary
    dirname = os.path.abspath( outputdir )
    if not os.path.isdir( dirname ):
        try:
            os.makedirs( dirname )
        except:
            raise ValueError("Cannot find or create directory '{}'".format( dirname ))

    # execute using qrsh
    if use_submit:
        qscript = "qsub -b y"
        stdoutname = "/dev/null"
        stderrname = "/dev/null"
        cmdelts = cmd.split(' ')
        i = 0
        while i < len(cmdelts):
            elt = cmdelts[i]
            if elt.startswith(">"):
                stdoutname = elt.lstrip('>')
                cmdelts[i] = ""
                if stdoutname == "":
                    i += 1
                    stdoutname = cmdelts[i]
                    cmdelts[i] = ""
            if elt.startswith("2>"):
                stderrname = elt[1:].lstrip('>')
                cmdelts[i] = ""
                if stderrname == "":
                    i += 1
                    stderrname = cmdelts[i]
                    cmdelts[i] = ""
            i += 1
        cmd = " ".join(cmdelts)
    else:
        stdoutname = "{}/qrsh.{}.stdout".format(dirname, hex(hash(cmd)))
        stderrname = "{}/qrsh.{}.stderr".format(dirname, hex(hash(cmd)))
        qscript = "qrsh"
    command = "{} -j y -cwd -now n -V -N {} -o {} -e {} {} {}".format(qscript, name, stdoutname, stderrname, qsub_config, cmd)
    if use_popen:
        popen = subprocess.Popen( command, shell=shell )
        return popen

    exitcode = subprocess.check_call( command, shell=shell )

    # clean up if successful
    if exitcode == 0 and not use_submit:
        try:
            os.unlink(stdoutname)
            os.unlink(stderrname)
        except OSError:
            pass

    # done
    return exitcode


def Popen(cmd, shell=True, outputdir=None, name=None, submit_local=False):
    return check_call(cmd, shell, outputdir, use_popen=True, submit_local=submit_local)


def submit(cmd, shell=True, outputdir=None, name=None, submit_local=False):
    return check_call(cmd, shell, outputdir, use_submit=True, submit_local=submit_local)
