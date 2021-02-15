# msn-model

Computational model of a striatial medium spiny neuron (MSN).

The model is implemented in [NEURON]. The associated Python files use
NEURON to build a MSN and simulate various experimental conditions. The
Python examples here were designed with the intention to provide some
understanding of the various steps involved in building and manipulating
a neuron model. After studying each example it should be possible to
combine these to simulate experiments and assess specific hypotheses.

2021 [Antonio González](mailto:antgon@cantab.net)

## Credits

The original NEURON model (msn/mechanisms/\*.mod), cell morphologies
(msn/morphologies/\*.swc), and best-fit parameters
(msn/parameters/\*.pkl) are from Lindroos and Hellgren Kotaleski
[Lindroos2020], available from ModelDB ([model ID
266775](http://modeldb.yale.edu/266775)). See also [Lindroos2018]. Some
of the Python code used to simulate MSNs is derived from code by the
same authors available [on
GitHub](https://github.com/robban80/striatal_SPN_lib); where applicable,
the documentation in the individual files acknowledges this fact.

## Requirements

* [Python 3 (>= 3.7)](https://www.python.org/)
* [Numpy](http://www.numpy.org/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [NEURON]

## How to use it

1. Install NEURON with Python; refer to [NEURON]'s website for details.
2. Run `nrnivmodl` inside the directory `msn/mechanisms` to compile the .mod
   files; see [NEURON]'s website.
3. Run the example scripts provided (e.g. `python example_1_build.py`).
   These files should be self explanatory.

## Links

* The [source code](https://github.com/antgon/msn-model) is on GitHub.

[NEURON]: https://neuron.yale.edu/neuron/

## References

[Lindroos2018]: Lindroos R, Dorst MC, Du K, Filipović M, Keller D,
Ketzef M, Kozlov AK, Kumar A, Lindahl M, Nair AG, Pérez-Fernández J,
Grillner S, Silberberg G & Hellgren Kotaleski J (2018). Basal ganglia
neuromodulation over multiple temporal and structural scales-simulations
of direct pathway MSNs investigate the fast onset of dopaminergic
effects and predict the role of Kv4.2. Front Neural Circuits 12, 3
(https://doi.org/10.3389/fncir.2018.00003).

[Lindroos2020]: Lindroos R & Hellgren Kotaleski J (2020). Predicting
complex spikes in striatal projection neurons of the direct pathway
following neuromodulation by acetylcholine and dopamine. Eur J Neurosci
(https://doi.org/10.1111/ejn.14891).
