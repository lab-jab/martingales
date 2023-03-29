## Replication files

### data folder  
- FX_rates.xlsx (exchange rates from the Federal Reserve of New York)
- The Wright (2001) data are included in the vrtest package

### code folder
- TS_gen.R (functions for generating returns)
- simulation.R (Monte Carlo simulations)
- FX.R (Application to exchange rates)

### convnets folder
Trained convnets in .h5 format. To load a model use: `m = load_model_hdf5("path_to_model/model_name.h5")`
- Wright2001.h5
- FRNY_weekly.h5
- FRNY_monthly.h5



