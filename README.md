# popdisp

**pop**ulation **disp**ersals: method to estimate allele frequencies in a region considering geographical locations of sampling sites, the nonequal number of samples in locations, and, most crucially, possible ways of dispersals within a region. It utilizes hierarchical Bayesian Modeling, compositional nature of frequencies and the Wright-Fisher model of drift.

## Initial data

Initial data contais two parts:  
(1) covariance matrix between samples based on the hypothetical dispersals within a region: `data/cov_mx`  
(2) allele frequencies in locations within regions `data\samples`

## Pipeline (running the test)

To demonstrate the test, please run:
* `estim_routes.py` to get estimates of allele frequencies
* `get_stat.py` to get statisticks for mcmc


## Requirements

To run **popdisp** methods, you need Python 3.4 or later. A list of required Python packages that the **popdisp** depends on, are in `requirements.txt`.  


## Authors

Anna Igolkina developed the **popdisp** package, [e-mail](mailto:igolkinaanna11@gmail.com).    


## License information
The **popdisp** package is open-sourced software licensed under the [MIT license](https://opensource.org/licenses/MIT).




