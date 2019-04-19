Welcome to Fish_N_CHIP

# What is Fish_N_CHIP ?

Fish_N_CHIP is a gathering of Snakemake workflows.

We use:
* python to parse command lines and launch snakfiles
* wrappers to support installation, testing, and execution

Nothing more.

# Requirements and installation

Each Snakemake workflow has its own requirements, please look at the Readme-s. However, you should be ok with:

1. Installation of [Conda](https://docs.anaconda.com/anaconda/install/linux/) for virtual environment managing

Download:`wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh`

Installation: `bash Anaconda3-2018.12-Linux-x86_64.sh`

2. Installation of [Datrie since it faile on python 3.7](https://github.com/pytries/datrie/issues/68), and Python 3.7 itself:

`conda create -n Python_3_7 datrie python=3.7`

4. Installation of snakemake and other dependencies

`python3.7 -m pip install snakemake pyyaml pyjson typing`

# Why Fish_N_CHIP ?

For two very important reasons:

* We use [Salmon](https://salmon.readthedocs.io/en/latest/), and [Feed the bears](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/) while working on RNA-Seq
* We have a lot of CHIP-Seq data in our labs
* We like Fish and Chips

Join us!
