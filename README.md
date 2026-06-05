# Semi-Leptonic Observable Predictions (SLOP)

[![DOI](https://zenodo.org/badge/1132340749.svg)](https://doi.org/10.5281/zenodo.18222257)

Repository for $b\to c\ell\nu$ (and some $b\to u\ell\nu$) observable predictions with varying form factors and Wilson coefficients. To be used to generate bands from fit results and to generate predictions for fits/fit models (e.g. for an unfolded fit).

You can find the code documentation online at [https://dvicoben.github.io/slophep/](https://dvicoben.github.io/slophep/). In case of questions, issues, or a particular request please open an issue on GitHub.

Predictions are made using `flavio` to compute the amplitudes and observables. ***Note that angular conventions may differ***. ***PDF/observable normalizations can also differ*** - literature and prediction tools vary in what factors are absorbed by the FFs, amplitudes, observables, decay rate and BR. For consistency it is best to look at rate-normalised observables.

# Reference
There is a zenodo doi for the repository. To cite SLOP please obtain your preferred citation from the [zenodo page](https://doi.org/10.5281/zenodo.18222257).

It is stressed that SLOP is a package that relies on or re-implements utilities from elsewhere. As such, any use of SLOP should also cite [flavio](https://flav-io.github.io/), [hammer](https://hammer.physics.lbl.gov/) and [eos](https://eoshep.org/). Please ensure you also cite those tools according to whatever is specified in their documentation.

Any usage of EvtGen correction weightings should cite the work of Florian Herren and Raynette van Tonder, [arXiv:2602.18378](https://arxiv.org/abs/2602.18378), and it is strongly encouraged to use the weightings they developed instead of SLOP if they serve your purpose: https://github.com/Herren/MCAmbulance.


# Requirements
Requirements are listed in `requirements.txt`
- Predictions use `flavio` to go from FFs and WCs to amplitudes and observables
- Internally also uses standard libraries (`numpy`, `matplotlib`), and `iminuit` for the (currently very limited) fitting functionality

# Set-up
## Quick
Ensure you are in a python environment with all requirements in `requirements.txt`, then
```
git clone https://github.com/dvicoben/slophep.git
cd slophep
source ./setup.sh
```
The script `setup.sh` simply appends `src/` to the `PYTHONPATH` so that contents therein will be found when running scripts. You will need to `source ./setup.sh` whenever you start a new terminal session.

## Using pip
In the python environment of your choice, 
```
git clone https://github.com/dvicoben/slophep.git
cd slophep
pip install -e .
```
which should install the package (`slophep`) and the required dependencies.


# Usage

- Some example scripts are in the `python` directory.
- Generation of predictions with varying FFs and WCs is shown in the minimal example `python/example_simple.py`. You can find a comparison of FF schemes in `python/compare_FFschemes.py`.
- Examples for generating error bands can be found in `python/example_fluctuations_obs.py` and `python/example_fluctuations_BR.py`.
- Additional FF schemes can be implemented - they need to inherit from the base `FormFactor` class and implement the `get_ff(q2)` method. See the [example](https://dvicoben.github.io/slophep/advanced.html#making-a-ff-scheme) in the documentation.
- There is a small example for comparisons with binned experimental results of the observables (`python/example_FLtau_comparison.py`)

There is some very limited fitting functionality that is in development. It should not be regarded as stable. For usage in fits, it is strongly encouraged that you do not rely on SLOP's fitting utilities at the moment. You may take the predictions and use them to build a fit model as best works for you.


# About the Predictions
Information about the predictions, as well as available form-factor schemes and decay modes, can be found in the online documentation [here](https://dvicoben.github.io/slophep/form-factors.html), or in its source [here](https://github.com/dvicoben/slophep/blob/master/docs/source/form-factors.md).

Note that there can be caveats for different decay modes and form-factors and it is strongly encouraged that users ensure that predictions they are using/plotting work as intended as I cannot cross-check/cover every use-case.


# TO DO
### Priority:
- [ ] Work on `PyPI` release
- [ ] Homogenise nomenclature of FF parameters for parameterisations with polynomial expansions
- [ ] Add ability to get $\langle J_i \rangle$ for a given binning scheme (as in the PDF methods) rather than need to get each individual bin
- [ ] Maybe move FF param defaults to some `.json` files? In particular for HPQCD this is a lot of parameters - largely a cosmetic thing and would like to keep everything readable from the class so maybe not

### Others
- [ ] Add some `cite` attirbute to return bib entries for each FF scheme - make bookkeeping easier for end-user
- [ ] Consider moving error handling to `gvar` rather than sampling of Gaussian? 
    - Moving to `gvar` could be problematic for non-gaussian errors - would need to consider how to deal with this.
    - In some tests for $B_s \to D_s^*$, obtained larger contours for the low $q^2$ in form factors compared to [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2) and what is obtained with LOAD_FIT.py. 
    - Even sampling directly from `gvar` (loading the `pydat` and using `gvar.sample`) produces errors different to those resulting from `gvar` arithmetic
    - The issues does not seem to be computation of the FFs since central values seem fine - need to cross-check directly by calculating similar fluctuations using functions in `LOAD_FIT.py`.
- [ ] Add ability to switch angular conventions
