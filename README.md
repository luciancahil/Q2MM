# Q2MM


## Conda

To install the conda env, run the following:

```
conda env create -f environment.yml
conda activate q2mm-openmm
```

## Testing Bayes optimization

Run the following:

```
bash run_bo.sh square_sin
```


## Testing the molecular mechanics

Run the following:

```
python get_graph.py "CC(C)=Cc1ccccc1c2ccccc2C=O"
bash run_bo.sh q2mm "CC(C)=Cc1ccccc1c2ccccc2C=O"
```