### Data analysis and visualisation pipeline for CosMx single-cell spatial transcriptomics (scST) readout

Version 0.9, M Chen, Imperial College, March 2024

### ModVoyager: 
scST analysis and visualisation
### DeSeq2Test: 
DEG analysis 

### Stopover: python implementation

Installation Guide (STopover)

## Installation and running
### 1. Python
#### Install conda environment and add jupyter kernel
```Plain Text  
  conda create -n STopover python=3.8
  conda activate STopover
  pip install git+https://github.com/bsungwoo/STopover.git
  pip install jupyter
  python -m ipykernel install --user --name STopover --display-name STopover
```
##### To annotate cells in image-based ST using TACCO:
```Plain Text  
  conda activate STopover
  pip install tacco
```
#### Run GUI for STopover (PyQt)
```Plain Text
  conda activate STopover
  python
  from STopover import app
  app.main()
```
