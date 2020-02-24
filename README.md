# DNA Re-replication model
 
A stochastic hybrid model for DNA re-replication, tailored for the complete fission yeast genome


## Prerequisites

* MATLAB R2011 or newer


## How to run model simulations

To perform a single DNA re-replication simulation, call `rereplicationalg` with the necessary impute arguments `copies`(genome amplification parameter C), `redistr` (binary variable indicating which model variant to use (1: LF, 0: UF)).
The example below simulates the LF model variant until an amplification level of 4C is reached:
```
copies=4;
redistr=1;
[Tfire,TPR,TSR,TSL,OS,evolution,lambdacurrent]=rereplicationalg(copies,redistr);
```

## License

This project is licensed under the terms of the MIT license.
