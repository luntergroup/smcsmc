from pytest import raises
from subprocess import run, CalledProcessError


def test_help():
    command = "smc2 --help"
    result = run(command, shell=True, check=True)
    assert result.returncode == 0


def test_version():
    """This should fail right now on CI as we have not installed smcsmc yet."""
    command = "smc2 -smcsmcpath build/smcsmc --version"
    result = run(command, shell=True, check=True, capture_output=True)
    assert result.returncode == 0
