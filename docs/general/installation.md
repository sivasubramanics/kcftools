# Installtation Guide

This guide explains the various ways to install `kcftools` on your system, including the recommended method via Bioconda, as well as options for downloading binaries or building from source. Choose the method that best fits your workflow and environment.

## 1. Install with Bioconda (Recommended)

The easiest and recommended way to install `kcftools` is through [Bioconda](https://anaconda.org/bioconda/kcftools). Make sure you have conda installed and set up with the Bioconda channel. 

If you don't have conda installed, you can follow the [Conda installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Or simply install Miniconda as follows:

    # create the directory for Miniconda (make sure to adjust the path if needed)
    $ mkdir -p ~/miniconda3

    # download and install Miniconda
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh

    # install Miniconda
    $ bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    # clean up the installer
    $ rm ~/miniconda3/miniconda.sh

    # After installing Miniconda, you need to initialize it. Run the following commands:
    $ source ~/miniconda3/bin/activate
    $ conda init --all

    # add Bioconda channel to conda configuration
    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge 

Once the conda is set up, you can install `kcftools` by running the following command:

    # it is optional to create a new conda environment for kcftools (but recommended)
    $ conda create -n kcftools -c bioconda kcftools
    $ conda activate kcftools
---

## 2. Install from Release Binary
If you prefer not to use conda, you can download the precompiled binary from the [KCFTOOLS releases page](https://github.com/sivasubramanics/kcftools/releases/tag/v0.1.0). Further follow the below steps to download and install:

    # download the latest release binary
    $ wget https://github.com/sivasubramanics/kcftools/releases/download/v0.1.0/kcftools-v0.1.0.zip
    # unzip the downloaded file
    $ unzip kcftools-v0.1.0.zip

    # move the binary to a directory in your PATH (optional)
    $ mv kcftools-v0.1.0/kcftools-0.1.0.jar $PATH/kcftools.jar
    # clean up the downloaded files
    $ rm kcftools-v0.1.0.zip

    # add HEAPSIZE variable and alias to the .bashrc file (optional)
    $ echo "export KCFTOOLS_HEAP_SIZE=4G" >> ~/.bashrc
    $ echo "alias kcftools='java -Xmx\$HEAPSIZE -jar \$PATH/kcftools.jar'" >> ~/.bashrc
    $ source ~/.bashrc
---

## 3. Install from Source

If you prefer to build from source, you can clone the repository and build it using Maven. This requires Java 11 or higher and Maven installed on your system.

    # clone the repository
    $ git clone https://github.com/sivasubramanics/kcftools
    $ cd kcftools

    # build the project using Maven
    $ mvn clean package

    # the built jar file will be in the target directory
    $ mv target/kcftools-0.1.0.jar $PATH/kcftools.jar
    
    # clean up the cloned repository (optional)
    $ cd ..
    $ rm -rf kcftools

    # add HEAPSIZE variable and alias to the .bashrc file (optional)
    $ echo "export KCFTOOLS_HEAP_SIZE=12G" >> ~/.bashrc
    $ echo "alias kcftools='java -Xmx\$HEAPSIZE -jar \$PATH/kcftools.jar'" >> ~/.bashrc
    $ source ~/.bashrc
---

## 4. Verify Installation

After installation, you can verify that `kcftools` is installed correctly by running:

    $ kcftools --version

This should display the version of `kcftools` you have installed.

---