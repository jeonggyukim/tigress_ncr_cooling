`tigress_ncr_cooling` is a stand-alone C++ program for calculating equilibrium photochemistry of the interstellar medium based on [Kim, Gong, Kim, & Ostriker 2023](https://ui.adsabs.harvard.edu/abs/2023ApJS..264...10K/abstract). The source code builds upon the photochemistry module implemented in the Athena-TIGRESS code (C version; to be updated) and infrastructure in [Athena++](https://github.com/PrincetonUniversity/athena).

To compile and run simulation:

```
./configure.py
make clean
make
cd inputs
python run_simulations.py
```

You may need to manually adjust input parameters (e.g., FUV, CR rate, metallicity, dust abundance) in `run_simulations.py`. Also pass `--help` after `configure.py` and 'run_simulations.py` to see other options.
