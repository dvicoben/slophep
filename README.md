# bd2dstlnu_angular

Repository for $B^0 \to D^*\mu\nu$ angular analysis predicitons.

To be used to generate bands from fit results and to generate predictions for fits/fit models (e.g. for an unfolded fit).


# Requirements

- Predictions are made using `flavio` to compute the amplitudes and observables - note that angular convention may differ
- Internally also uses standard libraries (`numpy`, `matplotlib`)

# Usage

- Some example scripts are in the `python` directory
- Additional FF schemes can be implemented - they need to inherit from `FormFactor` and implement the `get_ff(q2)` method, returning FFs in the lattice convention $V, A_0, A_1, A_{12}, T_1, T_2, T_{23}$. See existing schemes for examples.

# TO DO

- [x] Clean-up FF implementations 
    - There are a lot of values in these FF classes that are hidden underneath the `get_ff(q2)` computation - plan to change this to be grouped under some class attribute
- [x] Make documentation
- [ ] Add bibliography 
    - Mainly for easy reference of implemented schemes, make it easier to cross-check/verify
    - Also add it to documentation
- [ ] Sampling 
    - Ideally would provide some covariance matrix for FFs (and WCs) to resample from and generate confidence bands
- [ ] Easier binning predictions 
    - At moment need to provide boundaries for single bin to get prediction for that bin, plan to improve this to return predictions for all the range from a provided binning scheme
- [ ] Add 1D projections