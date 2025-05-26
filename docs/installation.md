## System Requirement

UMIche is advised to be installed within a conda environment, which makes it easier to work on multiple platforms, such as :material-microsoft-windows: Windows (partial), :simple-apple: Mac, and :material-linux: Linux. Owing to the exclusivity of Pysam to the Linux and Mac environments, UMIche does not work with BAM-related analysis in the Windows system.

!!! info "Note"

    Please note that starting from version `0.1.5`, **pysam** is no longer a required dependency when installing umiche. This means you can install umiche on any operating system. If you need to use **pysam**, you should install umiche on a non-Windows system and then install **pysam** separately. Please refer to [the documentation](https://pysam.readthedocs.io/en/latest/installation.html) :fontawesome-solid-arrow-up-right-from-square: for instructions on how to install it.

## :simple-pypi: PyPI (**highly recommended**, see [the latest version](./changelog.md))

[umiche homepage on PyPI](https://pypi.org/project/umiche)

!!! info "Note"

    Please make sure to use the latest version of umiche, as earlier versions may contain bugs. If you do not include the `--upgrade` flag during installation, you might encounter issues.


```shell
# create a conda environment
conda create --name umiche python=3.11

# activate the conda environment
conda activate umiche

# the latest version
pip install umiche --upgrade
```

## :simple-anaconda: Conda

[umiche homepage on Anaconda](https://anaconda.org/Jianfeng_Sun/umiche)

```shell
# create a conda environment
conda create --name umiche python=3.11

# activate the conda environment
conda activate umiche

# the latest version
conda install jianfeng_sun::umiche
```


## :fontawesome-brands-docker: Docker

[umiche homepage on Docker](https://hub.docker.com/r/2003100127/umiche)

```shell
docker pull 2003100127/umiche
```


## :octicons-mark-github-16: Github

[umiche homepage on Github](https://github.com/2003100127/umiche)

```shell
# create a conda environment
conda create --name umiche python=3.11

# activate the conda environment
conda activate umiche

# create a folder
mkdir project

# go to the folder
cd project

# fetch UMIche repository with the latest version
git clone https://github.com/2003100127/umiche.git

# enter this repository
cd umiche

# do the following command
pip install .
# or
python setup.py install
```
