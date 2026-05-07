# Using SLOP

## Getting a Prediction

To get a prediction, you must initialise a prediction object for the decay mode of interest and specify the form-factor scheme. You can then calculate a variety of observables. You can check available observables for available decay modes in the Python API.

```{code-block} python

from slophep.Predictions.Observables import BdToDstEllNuPrediction
from slophep.Predictions.FormFactorsBToV import BdToDstFF

pred = BdToDstEllNuPrediction("mu", "mu", BdToDstFF.CLN)
angobs = pred.J(q2=5.0)
print(angobs)
```

Where the example above computes rate-normalised angular observables in $B^0 \to D^{*}\mu\nu$ at $q^2 = 5.0$ $\text{GeV}^2$, using CLN form-factors, and prints out:

```
{'1s': 0.3805129434370365, '1c': 0.4933741604586339, '2s': 0.12646033489965114, '2c': -0.48972052780118125, '6s': -0.3912727821568956, '6c': 0.0035799980114745145, 3: -0.16082906069974, 4: 0.3182936823105255, 5: -0.3021657916707698, 7: -0.0, 8: 0.0, 9: -0.0}
```

You may change the form-factor parameters using the `pred.set_ff` and providing a dictionary with different parameters:

```{code-block} python

hqet2 = {
    "RhoSq" : 1.122,
    "h_A1"  : 0.908,
    "R1"    : 1.270,
    "R2"    : 0.852,
    "R0"    : 1.15
}
pred.set_ff(hqet2)
```

Similarly, you can specify Wilson coefficients using `obs.set_wc`:

```{code-block} python

wcoeffs = {
    "CVL_bcmunumu" : 0.0, 
    "CVR_bcmunumu" : 0.5,
    "CSL_bcmunumu" : 0.0,
    "CSR_bcmunumu" : 0.0,
    "CT_bcmunumu"  : 0.0
}
pred.set_wc(wcoeffs)
```

These are by default in the [WET `flavio` basis](https://wcxf.github.io/bases.html), but you may specify a different basis 
by through the additional `eft` and `basis` arguments, e.g. `pred.set_wc(wcoeffs, eft="WET", basis="flavio")`.


## Producing Errorbands

