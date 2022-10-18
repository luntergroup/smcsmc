from importlib.metadata import entry_points
from setuptools import setup

setup(
    name="smcsmc",
    version="1.1.0",
    description="Demographic Inference using Particle Filters",
    url="https://github.com/luntergroup/smcsmc",
    author="Chris Cole, Donna Henderson, Sha (Joe) Zhu, Gerton Lunter",
    author_email="ccole@well.ox.ac.uk, donnah@well.ox.ac.uk, sha.joe.zhu@gmail.com, gerton.lunter@imm.ox.ac.uk",
    license="CC-BY-4.0",
    packages=["smcsmc"],
    entry_points={
        "console_scripts": [
            "smc2 = smcsmc.cli:smcsmc_main",
        ]
    },
    install_requires=["numpy", "pandas", "msprime", "tqdm"],
    zip_safe=False,
)
