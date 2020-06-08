# DNA re-replication model
 
This repository contains the source code of our [preprint](https://www.biorxiv.org/content/10.1101/2020.03.30.016576v1) ***In silico* analysis of DNA re-replication across a complete genome reveals cell-to-cell heterogeneity and genome plasticity**. A complete mathematical description of the model can be also found [here](../master/DNA_rereplication_model.pdf). If you find our code useful please consider citing:

``` 
@article{rapsomaniki_silico_2020,
title = {In silico analysis of {DNA} re-replication across a complete genome reveals cell-to-cell 
heterogeneity and genome plasticity},
doi = {10.1101/2020.03.30.016576},
journal = {bioRxiv},
author = {Rapsomaniki, Maria Anna and Maxouri, Stella and Garrastacho, Manuel Ramirez and Nathanailidou, 
Patroula and Giakoumakis, Nickolaos Nikiforos and Taraviras, Stavros and Lygeros, John and Lygerou, Zoi},
month = mar,
year = {2020},
}
```


## Prerequisites

* MATLAB R2011 or newer


## How to run model simulations

To perform a single DNA re-replication simulation, call `rereplicationalg` with the necessary impute arguments:
- `copies`: genome amplification parameter *C*, 
- `redistr`: binary variable indicating if the *LF* (1) or the *UF* (0) model variant is simulated.

The example below simulates the LF model variant until an amplification level of `4C` is reached:
```
copies=4;
redistr=1;
[Tfire,TPR,TSR,TSL,OS,evolution,lambdacurrent]=rereplicationalg(copies,redistr);
```

## How to reproduce the preprint figures
In the MATLAB scripts `plots_paper.m` and `plots_supplementary.m` you can find the code needed to recreate all figures from our preprint and supplementary material. Please make sure to modify the path (first line of each cell) to match your local path where the simulation results were stored. 

## License

This project is free software, licensed under the terms of the MIT license.
