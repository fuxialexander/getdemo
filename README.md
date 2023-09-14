---
title: GET
app_file: app/main.py
colorFrom: gray
colorTo: red
sdk: docker
license: cc-by-nc-4.0
pinned: false
---


# Data preparation
Put the data in the following structure in the root directory of the project.
```bash
data
├── sequences
│   └── causal
│       ├── MECP2_TFAP2A
│       ├── PRDM1_SMAD2
│       └── TAF1_ZFX
└── structures
    ├── causal
    │   ├── MECP2_TFAP2A
    │   ├── PRDM1_SMAD2
    │   └── TAF1_ZFX
    └── homodimer
        ├── PRDM1
        ├── SMAD2
        ├── TAF1
        └── ZFX
```

# Installation
```bash
git clone --recursive git@github.com:fuxialexander/getdemo.git
cd getdemo
docker pull fuxialexander/getdemo:latest
docker run -it -v "/path/to/data:/data" --rm -p 7681:7681 fuxialexander/getdemo
or
singularity run  -w --bind /manitou/pmg/users/xf2217/getdemo:/app --bind /manitou/pmg/users/xf2217/demo_data:/data --bind /pmglocal/xf2217/tmp:/tmp --no-home --pwd /app getdemo
```
The gradio interface will be available at http://127.0.0.1:7681, a sharable link will be printed in the terminal.

# Build
```bash
git clone --recursive git@github.com:fuxialexander/getdemo.git
cd getdemo
docker build -t getdemo .
docker run -it -v "/path/to/data:/data" --rm -p 7681:7681 getdemo
```
