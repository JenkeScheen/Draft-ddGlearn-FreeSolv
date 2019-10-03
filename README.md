# Draft repository for learning ddGoffset values for FreeSolv molecules

**pFP/**: script to compute atom-pair fingerprints per datapoint.

**dPLEC/**: script to compute close contact fingerprints between ligands and proteins (irrelevant for freesolv).

**dFEAT/**: script to compute molecular properties per datapoint.

**datasets/**: folder containing all datasets to train on. Also contains the python script that compiles them.

**SVM/**: contains the script to run the training protocol. Might have to remove tensorflow imports.





**General workflow**: compile all molecules into **datasets/input/**, then run the python scripts for all three feature generation types. Once these have finished running, run compile_datasets.py in **datasets/** to combine all feature sets into different combinations. Once the datasets have been compiled, the builder script in **SVM/** can be run to start training; the other script plots convergence (this needs to be rewritten as it assumes splits per congeneric series). 

