`tigress_ncr_cooling` is a stand-alone C++ program for calculating equilibrium photochemistry of the interstellar medium based on [Kim, Gong, Kim, & Ostriker (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJS..264...10K/abstract). The source code builds upon the photochemistry module implemented in the Athena-TIGRESS code (written C; to be rewritten in C++) and infrastructure in [Athena++](https://github.com/PrincetonUniversity/athena).

To compile the code and run simulations:

```
./configure.py
make clean
make
cd inputs
python run_simulations.py
```

You may need to manually adjust input parameters (e.g., FUV, CR rate, metallicity, dust abundance) in `run_simulations.py`. Also pass `--help` to python scripts to see other options.

See [notebook](notebook) for results of some example calculations.
