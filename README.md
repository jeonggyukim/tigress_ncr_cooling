tigress_ncr_cooling is a stand-alone c++ program for calculating equilibrium photochemistry of the interstellar medium based on J.-G. Kim et al. (2023). The source code builds upon the photochemistry module implemented in the Athena-TIGRESS code and infrastructure in [Athena++](https://github.com/PrincetonUniversity/athena).

To compile and run simulation:

```
./configure.py
make clean
make
cd inputs
python run_simulations.py
```

The output files will be generated in bin directory.