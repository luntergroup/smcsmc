# Running through Docker

Docker provides containerization to enable applications to be deployed on cloud compute services. We have provided a recipe to build a simple Docker container for SMCSMC which can be either used on its own to launch custom jobs on AWS or GCP or as the starting point for more complex workflow configurations. 

:::{warning}
This is an advanced topic that requires knowledge of Docker and cloud compute services. To use SMCSMC with Docker, you must first install Docker. Follow the [official documentation](https://docs.docker.com/get-docker/) on instructions on how to do so. 
:::

This recipe can be found in the `./Dockerfile` of the [main repository](https://github.com/luntergroup/smcsmc) for SMCSMC. We host this docker image on Docker Hub with the following reserved tags along with release-specific tags:

| Description           | Tag           |
|-----------------------|---------------|
| Latest stable version | `latest`      |
| Development           | `test`        |
| Release-specific      | e.g. `v1.0.x` |

## Running commands within a Docker image

Simple commands can be piped directly to the docker image.  This will pull the docker image if it is not already installed locally. 

```
docker run chris1221/smcsmc:latest smc2 -h
```

More complex commands, as are typical when analysing real data, require either

1. A workflow manager such as Cromwell, Snakemake, or NextFlow to automate the calling of specific SMCSMC stages or
2. Custom configurations for your particular use case. An example of one is given below.


:::{note}
This example recreates one of the test cases in `test/test_conversions.py` within Docker 
:::


As an example, suppose we want to convert a VCF file to a Seg file (required for input to SMCSMC) using this docker image. One way to do this would be to create a directory with:

* The VCF files 
* A Python script to run SMCSMC 

This directory could then be mounted and the python script could be run within the context of the docker interpreter, which has `smcsmc` already installed. 

As an example, from the root of the `smcsmc` Github repo:

```sh
git clone https://github.com/luntergroup/smcsmc.git smcsmc_git
cd smcsmc_git 
cp -r test/data/ mount_directory 
```

And create the following python script in `mount_directory/convert_vcf.py`

```py
import smcsmc 

test_vcf = "/app/mount_directory/test.vcf.gz"
test_mask = "/app/mount_directory/test.bed.gz"
chr = 2

sample_1 = "ID1"
sample_2 = "ID2"

smcsmc.vcf_to_seg(
    [(test_vcf, sample_1), (test_vcf, sample_2)],
    "/app/mount_directory/testNoMask.seg.gz",
    tmpdir="test/out/tmp",
    key="testNoMask",
    chroms=[chr],
)
```

:::{note}
We can send shell commands directly to the image without starting it first because the default `entrypoint` is set to be `/bin/bash`; if your image has a non-standard `entrypoint` you must override the entrypoint to be `/bin/bash`. 
:::

Now run the docker container with the command to invoke the Python script, mounting the `mount_directory` to `/app/mount_directory`

```sh
echo python /app/mount_directory/convert_vcf.py | docker run -i -v "$(pwd)"/mount_directory:/app/mount_directory smcsmcgit:latest 
```

Because mounting is a bidirectional sync with the docker container, the resulting VCF files are available in the `mount_directory` upon completion of the command. This illustrates a simple example of running code from SMCSMC from within a Docker image. 