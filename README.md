# bd2dstlnu_angular

Repository for $B^0 \to D^*\mu\nu$ angular analysis predicitons.

To be used to generate bands from fit results and to generate predictions for fits/fit models (e.g. for an unfolded fit).


# Requirements

- Predictions are made using `flavio` to compute the amplitudes and observables - note that angular convention may differ
- Internally also uses standard libraries (`numpy`, `matplotlib`)


# Usage

- Some example scripts are in the `python` directory
- Additional FF schemes can be implemented - they need to inherit from `FormFactor` and implement the `get_ff(q2)` method, returning FFs in the lattice convention $V, A_0, A_1, A_{12}, T_1, T_2, T_{23}$. See existing schemes for examples.

## Implemented FFs

| FF Scheme | Implementation | Notes | Refs. |
|-----------|----------------|-------|-------|
| CLN       | flavio based   | Tensor FFs are obtained as in flavio, using eqn. 11 in [arXiv:1503.05534](https://arxiv.org/abs/1503.05534). Use with care for BSM predictions. Defaults are set to HAMMER values. | [arXiv:hep-ph/9712417v1](https://arxiv.org/abs/hep-ph/9712417), [arXiv:1203.2654](https://arxiv.org/abs/1203.2654), [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/clnexp.py), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarCLN.cc?ref_type=tags)
| BGL       | Hammer based   | Hammer divides by an additional factor of $\eta_{EW}V_{cb}$. Translation to $V, A_i$ obtained from EOS. Tensor FFs $T_i = 0$, use with care for BSM predictions. Defaults are set from [arXiv:1707.09509](https://arxiv.org/abs/1707.09509). |  [arXiv:hep-ph/9705252](https://arxiv.org/abs/hep-ph/9705252), [arXiv:1707.09509](https://arxiv.org/abs/1707.09509), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgl1997-impl.hh)|
| BSZ       | flavio based   | Resonances used taken from [arXiv:1811.00983](https://arxiv.org/abs/1811.00983). Should match EOS implementation (`BSZ2015` which is default $B\to D^*$ FF scheme in EOS). Defaults are set to EOS values (see [EOS docs](https://eoshep.org/doc/reference/parameters.html#parameters-in-b-to-v-form-factor-parametrizations)). | [arXiv:1503.05534](https://arxiv.org/abs/1503.05534), [arXiv:1811.00983](https://arxiv.org/abs/1811.00983), [flavio source](https://github.com/flav-io/flavio/blob/master/flavio/physics/bdecays/formfactors/b_v/bsz.py), [eos source](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bsz2015-impl.hh) |
| BLPR      | Hammer based   | Correspondence to $V, A_i, T_i$ obtained from Appendix B in [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) / Eqns. 38-39 in [arXiv:1309.0301](https://arxiv.org/abs/1309.0301) and similar parametrisation in eos (see [EOS BGJvD implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh)). Defaults are set to HAMMER values. | [arXiv:1703.05330](https://arxiv.org/abs/1703.05330), [arXiv:1908.09398](https://arxiv.org/abs/1908.09398), [hammer source](https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBLPR.cc?ref_type=tags) |


# TO DO

- [x] Clean-up FF implementations 
    - There are a lot of values in these FF classes that are hidden underneath the `get_ff(q2)` computation - plan to change this to be grouped under some class attribute
- [x] Make documentation
- [ ] Implement futher FF schemes
    - BGJvD (see [eos implementation](https://github.com/eos/eos/blob/v1.0.13/eos/form-factors/parametric-bgjvd2019-impl.hh), [arXiv:1912.09335](https://arxiv.org/abs/1912.09335)) and BLPRXP (see [hammer implementation](https://gitlab.com/mpapucci/Hammer/-/blob/v1.4.1/src/FormFactors/BLPRXP/FFBtoDstarBLPRXPVar.cc), [arXiv:2206.11281](https://arxiv.org/abs/2206.11281)), are BLPR-like schemes with subsubleading contributions ($\mathcal{O}(\varepsilon_c^2)$, $\mathcal{O}(\varepsilon_b\varepsilon_c)$)
    - `flavio` default FF is also BLPR-like and likely corresponds to [arXiv:1908.09398](https://arxiv.org/abs/1908.09398) - may be similar to EOS's BGJvD
- [x] Add bibliography 
    - Mainly for easy reference of implemented schemes, make it easier to cross-check/verify
    - Also add it to FF implementation documentation/docstrings
- [ ] Sampling 
    - Ideally would provide some covariance matrix for FFs (and WCs) to resample from and generate confidence bands
- [ ] Easier binning predictions 
    - At moment need to provide boundaries for single bin to get prediction for that bin, plan to improve this to return predictions for all the range from a provided binning scheme
- [ ] Add 1D projections

