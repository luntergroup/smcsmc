"""
Utilities for scripts
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import functools
import os
import shlex
import subprocess
import sys
import time

import humanize
import requests
import yaml
import pysam


def log(message):
    print(message)


class Timed(object):
    """
    Decorator that times a method, reporting runtime at finish
    """
    def __call__(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self.start = time.time()
            result = func(*args, **kwargs)
            self.end = time.time()
            self._report()
            return result
        return wrapper

    def _report(self):
        delta = self.end - self.start
        timeString = humanize.time.naturaldelta(delta)
        log("Finished in {} ({} seconds)".format(timeString, delta))


def runCommandSplits(splits, silent=False):
    """
    Run a shell command given the command's parsed command line
    """
    if silent:
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(splits, stdout=devnull, stderr=devnull)
    else:
        subprocess.check_call(splits)


def runCommand(command, silent=False):
    """
    Run a shell command
    """
    splits = shlex.split(command)
    runCommandSplits(splits, silent=silent)


def getAuthValues(filePath='scripts/auth.yml'):
    """
    Return the script authentication file as a dictionary
    """
    return getYamlDocument(filePath)


def getYamlDocument(filePath):
    """
    Return a yaml file's contents as a dictionary
    """
    with open(filePath) as stream:
        doc = yaml.load(stream)
        return doc


