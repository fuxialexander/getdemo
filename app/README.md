# Installation
```bash
docker pull fuxialexander/getdemo:latest
docker run -it -v "/path/to/data:/data" --rm -p 7681:7681 fuxialexander/getdemo
```
The gradio interface will be available at http://127.0.0.1:7681, a sharable link will be printed in the terminal.

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