# bd2dstlnu_angular

Repository for $B^0 \to D^*\ell\nu$ angular analysis predicitons.

To be used to generate bands from fit results and to generate predictions for fits/fit models (e.g. for an unfolded fit).

Predictions are made using `flavio` to compute the amplitudes and observables. ***Note that angular conventions may differ***. ***PDF/observable normalizations can also differ*** - literature and prediction tools vary in what factors are absorbed by the FFs, amplitudes, observables, decay rate and BR. For consistency it is best to look at rate-normalised observables.

# Requirements
Requirements are listed in `requirements.txt`
- Predictions use `flavio` to go from FFs and WCs to amplitudes and observables
- Internally also uses standard libraries (`numpy`, `matplotlib`), and `iminuit` for the (currently very limited) fitting functionality

# Set-up
## Quick
Ensure you are in a python environment with all requirements in `requirements.txt`, then
```
git clone https://gitlab.cern.ch/dvicoben/bd2dstmunu.git
cd bd2dstlnu_angular
source ./setup.sh
```
The script `setup.sh` simply appends `src/` to the `PYTHONPATH` so that contents therein will be found when running scripts. You will need to `source ./setup.sh` whenever you start a new terminal session.

## Using pip
In the python environment of your choice, 
```
git clone https://gitlab.cern.ch/dvicoben/bd2dstmunu.git
cd bd2dstlnu_angular
pip install -e .
```
which should install the package (`b2dstlnu`) and the required dependencies.


# Usage

- Some example scripts are in the `python` directory.
- Generation of predictions with varying FFs and WCs is shown in the minimal example `python/example_simple.py`. You can find a comparison of FF schemes (with some example of plotting functionality) in `python/compare_FFschemes.py`.
- Examples for generating error bands can be found in `python/example_fluctuations_obs.py` and `python/example_fluctuations_BR.py`.
- Additional FF schemes can be implemented - they need to inherit from `FormFactor` and implement the `get_ff(q2)` method, returning FFs in the basis $V, A_0, A_1, A_{12}, T_1, T_2, T_{23}$. See existing schemes for examples.
- There are some preliminary scripts for fits (largely illustrative), `python/test_coef_fit.py` and `python/test_FF_fit.py`. Currently working on more optimised fitting interface.
- Access (available) documentation in `docs/build/html`, open `index.html` in your preferred browser.

## Implemented FF Schemes

| FF Scheme | Implementation | Notes | Refs. |
|-----------|----------------|-------|-------|
| CLN       | flavio based   | Tensor FFs are obtained as in flavio, using eqn. 11 in [arXiv:1503.05534](https://arxiv.org/abs/1503.05534). Use with care for BSM predictions. Defaults are set to HAMMER values. | [arXiv:hep-ph/9712417v1](https://arxiv.org/abs/hep-ph/9712417), [arXiv:1203.2654](https://arxiv.org/abs/1203.2654), [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarCLN.cc?ref_type=tags)
| BGL       | Hammer based   | Hammer divides by an additional factor of $\eta_{EW}V_{cb}$. Translation to $V, A_i$ obtained from EOS. Tensor FFs $T_i = 0$, use with care for BSM predictions. Defaults are set from [arXiv:1707.09509](https://arxiv.org/abs/1707.09509). |  [arXiv:hep-ph/9705252](https://arxiv.org/abs/hep-ph/9705252), [arXiv:1707.09509](https://arxiv.org/abs/1707.09509), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgl1997-impl.hh)|
| BSZ       | flavio based   | FF from fit to LCSR + zero recoil lattice. Resonances used taken from [arXiv:1811.00983](https://arxiv.org/abs/1811.00983). Should match EOS implementation (`BSZ2015` which is default $B\to D^*$ FF scheme in EOS). Defaults are set to EOS values (see [EOS docs](https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-v-form-factor-parametrizations)). | [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [arXiv:1811.00983](https://arxiv.org/abs/1811.00983), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/bsz.py), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bsz2015-impl.hh) |
| BLPR      | Hammer based   | Correspondence to $V, A_i, T_i$ obtained from Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) / Eqns. 38-39 in [arXiv:1309.0301](https://arxiv.org/abs/1309.0301) and similar parametrisation in eos (see [EOS BGJvD implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh)). Defaults are set to HAMMER values. | [arXiv:1703.05330](https://arxiv.org/abs/1703.05330), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBLPR.cc?ref_type=tags) |
| HPQCD     | From ancillary files in [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2) | FF from fit to non-zero recoil lattice QCD in [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2), as described in Sec. IV B. | [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2) |


# TO DO
### Priority:
- [ ] Homogenise nomenclature of FF parameters for parameterisations with polynomial expansions
- [ ] Add ability to get $\langle J_i \rangle$ for a given binning scheme (as in the PDF methods) rather than need to get each individual bin
- [ ] Maybe move FF param defaults to some `.json` files? In particular for HPQCD this is a lot of parameters - largely a cosmetic thing and would like to keep everything readable from the class so maybe not

### Others
- [ ] Add some `cite` attirbute to return bib entries for each FF scheme - make bookkeeping easier for end-user
- [ ] Fitting interface for $\langle J_i \rangle$ fits
- [ ] Fitting interface for FF fits
- [ ] For fitting: Optimise binned PDF predictions to avoid unnecessary re-calculations of angular integrals unless explicitly requested
- [ ] Add 1D projections
- [ ] Implement futher FF schemes
    - BGJvD (see [eos implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh), [arXiv:1912.09335](https://arxiv.org/abs/1912.09335)) and BLPRXP (see [hammer implementation](https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BLPRXP/FFBtoDstarBLPRXPVar.cc), [arXiv:2206.11281](https://arxiv.org/abs/2206.11281)), are BLPR-like schemes with subsubleading contributions ($\mathcal{O}(\varepsilon_c^2)$, $\mathcal{O}(\varepsilon_b\varepsilon_c)$)
    - `flavio` default FF is also BLPR-like and likely corresponds to [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) - may be similar to EOS's BGJvD
- [ ] Add ability to switch angular conventions

### Done
- [x] Clean-up FF implementations 
    - There are a lot of values in these FF classes that are hidden underneath the `get_ff(q2)` computation - plan to change this to be grouped under some class attribute
- [x] Make documentation
- [x] Add bibliography 
    - Mainly for easy reference of implemented schemes, make it easier to cross-check/verify
    - Also add it to FF implementation documentation/docstrings
- [x] Easier binning predictions
    - At moment need to provide boundaries for single bin to get prediction for that bin, plan to improve this to return predictions for all the range from a provided binning scheme
- [x] Errorbands / Sampling 
    - Ideally should be able provide some covariance matrix for FFs (and WCs) to resample from and generate confidence bands
- [x] Update repo to work like a python module
    - Need the `.toml` and requirements for an easy `pip` install
- [x] Add base methods in `FormFactor` to obtain FFs in different basis from currently mandatory $A_i$, i.e.
    - [x] Common HQET basis $h_{A_i}$
    - [x] Usual BGL FFs $g/f/F_1/F_2$ 
    - [x] Standard CLN $h/R_1/R_2/R_0$


