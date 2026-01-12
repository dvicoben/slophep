# Semi-Leptonic Observable Predictions (SLOP)

[![DOI](https://zenodo.org/badge/1132340749.svg)](https://doi.org/10.5281/zenodo.18222257)

Repository for $b\to c\ell\nu$ (and some $b\to u\ell\nu$) observable predictions with varying form factors and Wilson coefficients.

To be used to generate bands from fit results and to generate predictions for fits/fit models (e.g. for an unfolded fit).

Predictions are made using `flavio` to compute the amplitudes and observables. ***Note that angular conventions may differ***. ***PDF/observable normalizations can also differ*** - literature and prediction tools vary in what factors are absorbed by the FFs, amplitudes, observables, decay rate and BR. For consistency it is best to look at rate-normalised observables.

In case of questions, issues, or a particular request you can contact me by mattermost or email (david.vico.benet@cern.ch).

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
- Additional FF schemes can be implemented - they need to inherit from `FormFactorBToV` and implement the `get_ff(q2)` method, returning FFs in the basis $V, A_0, A_1, A_{12}, T_1, T_2, T_{23}$. See existing schemes for examples.
- There is a small example of how to handle comparisons with binned experimental results of the observables (`python/example_FLtau_comparison.py`)
- There are some preliminary scripts for fits (largely illustrative), `python/test_coef_fit.py` and `python/test_FF_fit.py`. Currently working on more optimised fitting interface.
- Access (available) documentation in `docs/build/html`, open `index.html` in your preferred browser.


# About the Predictions
For predicting absolute values (e.g. the BR of a decay mode over $q^2$), SLOP uses hadron masses, $G_F$, CKM matrix elements, and hadron lifetimes from `flavio`'s `default_parameters`. This can result in slight differences with respect to predictions from elsewhere (e.g. EOS) even if the FF and amplitude/observable computation is the same due to these differing factors. For rate normalised observables or ratios these factors should cancel and there should be no difference.

Also important to note: Wilson coefficients are defined at a particular scale. In SLOP the scale is set to $\mu = 4.8$ by default, but this can be changed when initialising the prediction. For plotting error bands this is particularly important as your fitter may not assume the same scale.

For $B \to P$ predictions, the decomposition of the decay rate, up to normalisation constant, follows the form:
$$\frac{\mathrm{d}\Gamma}{\mathrm{d}q^2 \mathrm{d}\cos\theta_\ell} \propto
a(q^2) + b(q^2)\cos\theta_\ell + c(q^2)\cos^2\theta_\ell
$$

For $B \to V$ predictions, the decomposition of the decay rate used throughout, up to a normalisation constant, follows the form:

$$\frac{\mathrm{d}\Gamma}{\mathrm{d}q^2 \mathrm{d}\cos\theta_V \mathrm{d}\cos\theta_\ell \mathrm{d}\chi} 
\propto \frac{9}{32\pi} \Bigg[ 
J_{1c}\cos^2\theta_V + J_{1s}\sin^2\theta_V
\\ +\left(J_{2c}\cos^2\theta_V + J_{2s}\sin^2\theta_V \right)\cos 2\theta_\ell
\\ +\left(J_{6c}\cos^2\theta_V + J_{6s}\cos^2\theta_V \right)\cos\theta_\ell
\\ +\left(J_3\cos 2\chi + J_9 \sin 2\chi \right)\sin^2\theta_\ell \sin^2\theta_V
\\ +\left(J_4 \cos\chi + J_8\sin\chi\right)\sin 2\theta_\ell \sin 2\theta_V
\\ + \left(J_5 \cos\chi + J_7 \sin\chi \right) \sin\theta_\ell \sin 2\theta_V
\Bigg]$$


# Available Decay Modes & FF Schemes
## $B \to \pi$

| FF Scheme | Notes | Refs. |
|-----------|-------|-------|
| BSZ       | Implementation reproduces flavio's. FF from fit to LCSR + zero recoil lattice. Resonances used taken from [arXiv:1811.00983](https://arxiv.org/abs/1811.00983). Should match EOS implementation (`BSZ2015`). Defaults are set to EOS values (see [EOS docs](https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-p-form-factor-parametrizations)). | [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [arXiv:1811.00983](https://arxiv.org/abs/1811.00983), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_p/bsz.py), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bsz2015-impl.hh) |

## $B \to D$

Hammer equivalents available: `CLN_Hammer`, `BGL_Hammer`, `BLPR_Hammer`. Note that the $B \to P$ basis in Hammer differs in $f_T$ from the one used in flavio and therefore in SLOP for the calculation of observables (by a factor of $m_B + m_P$). Because of this, take care when comparing $f_T$ (and any observables it affects) from `BLPR_Hammer` in $B\to P$ and other FF schemes (`BGL_Hammer` and `CLN_Hammer` are SM only so this is not an issue as $f_T = 0$).

For BGL parameterizations, references vary in the resonance masses used in the Blaschke factors as well as the $\chi$'s used in the outer functions. A one-to-one correspondence is not expected unless all of these match with SLOP, even if one sets the same expansion coefficients.

| FF Scheme | Notes | Refs. |
|-----------|-------|-------|
| CLN       | Implementation reproduces hammer's. Tensor FF $f_T(q^2) = 0$. Use with care for BSM predictions. Defaults are set to HAMMER values. | [arXiv:hep-ph/9712417v1](https://arxiv.org/abs/hep-ph/9712417), [arXiv:1203.2654](https://arxiv.org/abs/1203.2654), [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDCLN.cc?ref_type=tags)
| BGL       | Implementation reproduces Hammer's. Tensor FF $f_T(q^2) = 0$. Use with care for BSM predictions. Defaults are set to HAMMER values. |  [arXiv:hep-ph/9705252](https://arxiv.org/abs/hep-ph/9705252), [arXiv:1707.09509](https://arxiv.org/abs/1707.09509), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBGL.cc?ref_type=tags) |
| BSZ       | Implementation reproduces flavio's. FF from fit to LCSR + zero recoil lattice. Resonances used taken from [arXiv:1811.00983](https://arxiv.org/abs/1811.00983). Should match EOS implementation (`BSZ2015`). Defaults are set to EOS values (see [EOS docs](https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-p-form-factor-parametrizations)). | [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [arXiv:1811.00983](https://arxiv.org/abs/1811.00983), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_p/bsz.py), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bsz2015-impl.hh) |
| BLPR      | Implementation reproduces hammer's. Correspondence to from $h_i$ to $f_i$ should follow Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398). Defaults are set to HAMMER values. | [arXiv:1703.05330](https://arxiv.org/abs/1703.05330), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDBLPR.cc?ref_type=tags) |
| BLPRXP (checks pending) | Implementation intended to reproduce Hammer's. Correspondence to from $h_i$ to $f_i$ should follow Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398). Defaults are set to HAMMER values. | [arXiv:2206.11281](https://arxiv.org/abs/2206.11281), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BLPRXP/FFBtoDstarBLPRXP.cc) |

## $B \to D^*$

Hammer equivalents available: `CLN_Hammer`, `BGL_Hammer`, `BLPR_Hammer`. `BGL_Hammer` has the additional $1/\eta_{EW}V_{cb}$ that is present in Hammer but not in SLOP's `BGL`.

For BGL parameterizations, references vary in the resonance masses used in the Blaschke factors as well as the $\chi$'s used in the outer functions. A one-to-one correspondence is not expected unless all of these match with SLOP, even if one sets the same expansion coefficients.

| FF Scheme | Notes | Refs. |
|-----------|-------|-------|
| ISGW2     | Implementation meant to reproduce Hammer's (see [source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarISGW2.cc)). Hammer itself is based on the EvtGen model. Correspondence of hammer's basis to flavio/SLOP basis using [Hammer Manual](https://hammer.physics.lbl.gov/HammerManual.pdf) (Appendix D.2) and Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398). | [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarISGW2.cc) |
| CLN       | Implementation reproduces flavio's. Tensor FFs are zero. There is a similar `CLN2` implentation where tensor FFs are obtained as in flavio, using eqn. 11 in [arXiv:1503.05534](https://arxiv.org/abs/1503.05534). Use with care for BSM predictions. Defaults are set to HAMMER values. | [arXiv:hep-ph/9712417v1](https://arxiv.org/abs/hep-ph/9712417), [arXiv:1203.2654](https://arxiv.org/abs/1203.2654), [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarCLN.cc?ref_type=tags)
| BGL       | Implementation reproduces Hammer's. Hammer divides by an additional factor of $\eta_{EW}V_{cb}$. Translation to $V, A_i$ obtained from EOS. Tensor FFs $T_i = 0$, use with care for BSM predictions. Defaults are set to HAMMER values divided by $\eta_{EW}V_{cb}$. |  [arXiv:hep-ph/9705252](https://arxiv.org/abs/hep-ph/9705252), [arXiv:1707.09509](https://arxiv.org/abs/1707.09509), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgl1997-impl.hh)|
| BSZ       | Implementation reproduces flavio's. FF from fit to LCSR + zero recoil lattice. Resonances used taken from [arXiv:1811.00983](https://arxiv.org/abs/1811.00983). Should match EOS implementation (`BSZ2015` which is default $B\to D^*$ FF scheme in EOS). Defaults are set to EOS values (see [EOS docs](https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-v-form-factor-parametrizations)). | [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [arXiv:1811.00983](https://arxiv.org/abs/1811.00983), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/bsz.py), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bsz2015-impl.hh) |
| BLPR      | Implementation reproduces hammer's. Correspondence to $V, A_i, T_i$ obtained from Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) / Eqns. 38-39 in [arXiv:1309.0301](https://arxiv.org/abs/1309.0301) and similar parametrisation in eos (see [EOS BGJvD implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh)). Defaults are set to HAMMER values. | [arXiv:1703.05330](https://arxiv.org/abs/1703.05330), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBLPR.cc?ref_type=tags) |
| BLPRXP (checks pending) | Implementation intended to reproduce Hammer's. Correspondence to $V, A_i, T_i$ obtained from Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) / Eqns. 38-39 in [arXiv:1309.0301](https://arxiv.org/abs/1309.0301) and similar parametrisation in eos (see [EOS BGJvD implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh)). Defaults are set to HAMMER values. | [arXiv:2206.11281](https://arxiv.org/abs/2206.11281), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BLPRXP/FFBtoDstarBLPRXP.cc) |
| HPQCD     | Implementation from ancillary files in [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2). FF from fit to non-zero recoil lattice QCD in [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2), as described in Sec. IV B. | [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2) |


## $B_s \to D_s^*$
Implementation in most cases is the same as $B \to D^*$ but with the appropriate meson masses. Note that [arXiv:1801.10468](https://arxiv.org/pdf/1801.10468) uses a different angular decomposition for $D^*\to D\gamma$ which is not available in SLOP.

Hammer equivalents available: `CLN_Hammer`, `BGL_Hammer`, `BLPR_Hammer`.

For BGL parameterizations, references vary in the resonance masses used in the Blaschke factors as well as the $\chi$'s used in the outer functions. A one-to-one correspondence is not expected unless all of these match with SLOP, even if one sets the same expansion coefficients.

| FF Scheme | Notes | Refs. |
|-----------|-------|-------|
| ISGW2     | Implementation meant to reproduce Hammer's (see [source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarISGW2.cc)). Hammer itself is based on the EvtGen model. Correspondence of hammer's basis to flavio/SLOP basis using [Hammer Manual](https://hammer.physics.lbl.gov/HammerManual.pdf) (Appendix D.2) and Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398). | [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarISGW2.cc) |
| CLN       | Implementation reproduces flavio's. Tensor FFs are zero. There is a similar `CLN2` implentation where tensor FFs are obtained as in flavio, using eqn. 11 in [arXiv:1503.05534](https://arxiv.org/abs/1503.05534). Use with care for BSM predictions. Defaults are set to HAMMER values. | [arXiv:hep-ph/9712417v1](https://arxiv.org/abs/hep-ph/9712417), [arXiv:1203.2654](https://arxiv.org/abs/1203.2654), [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarCLN.cc?ref_type=tags)
| BGL       | Implementation reproduces Hammer's. Hammer divides by an additional factor of $\eta_{EW}V_{cb}$. Translation to $V, A_i$ obtained from EOS. Tensor FFs $T_i = 0$, use with care for BSM predictions. Defaults are set to HAMMER values divided by $\eta_{EW}V_{cb}$. |  [arXiv:hep-ph/9705252](https://arxiv.org/abs/hep-ph/9705252), [arXiv:1707.09509](https://arxiv.org/abs/1707.09509), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags) |
| BLPR      | Implementation reproduces hammer's. Correspondence to $V, A_i, T_i$ obtained from Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) / Eqns. 38-39 in [arXiv:1309.0301](https://arxiv.org/abs/1309.0301) and similar parametrisation in eos (see [EOS BGJvD implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh)). Defaults are set to HAMMER values. | [arXiv:1703.05330](https://arxiv.org/abs/1703.05330), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBLPR.cc?ref_type=tags) |
| BLPRXP (checks pending) | Implementation intended to reproduce Hammer's. Correspondence to $V, A_i, T_i$ obtained from Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) / Eqns. 38-39 in [arXiv:1309.0301](https://arxiv.org/abs/1309.0301) and similar parametrisation in eos (see [EOS BGJvD implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh)). Defaults are set to HAMMER values. | [arXiv:2206.11281](https://arxiv.org/abs/2206.11281), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BLPRXP/FFBtoDstarBLPRXP.cc) |
| HPQCD     | Implementation from ancillary files in [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2). FF from fit to non-zero recoil lattice QCD in [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2), as described in Sec. IV B. | [arXiv:2304.03137v2](https://arxiv.org/abs/2304.03137v2) |

## $\Lambda_b \to \Lambda_c$
Decay rate and observables follow definitions in [arXiv:1907.12554](https://arxiv.org/abs/1907.12554) as in EOS.

Note: SLOP by default uses $\alpha_{-}^{\Lambda_c} = -0.786$

| FF Scheme | Notes | Refs. |
|-----------|-------|-------|
| DKMR      | Implementation meant to reproduce EOS's (see [source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-dkmr2017-impl.hh)), which implements [arXiv:1702.02243](https://arxiv.org/abs/1702.02243) FFs and calculates observables as per [arXiv:1907.12554](https://arxiv.org/abs/1907.12554).  | [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-dkmr2017-impl.hh), [arXiv:1702.02243](https://arxiv.org/abs/1702.02243), [arXiv:1907.12554](https://arxiv.org/abs/1907.12554) |

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
- [ ] Implement more FF schemes
    - BGJvD (see [eos implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh), [arXiv:1912.09335](https://arxiv.org/abs/1912.09335)) and BLPRXP (see [hammer implementation](https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BLPRXP/FFBtoDstarBLPRXPVar.cc), [arXiv:2206.11281](https://arxiv.org/abs/2206.11281)), are BLPR-like schemes with subsubleading contributions ($\mathcal{O}(\varepsilon_c^2)$, $\mathcal{O}(\varepsilon_b\varepsilon_c)$)
    - `flavio` default FF is also BLPR-like and likely corresponds to [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) - may be similar to EOS's BGJvD
- [ ] Add ability to switch angular conventions
