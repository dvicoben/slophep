# Advanced Use

## Making a FF Scheme

It may be that there is a particular form-factor scheme you want to implement. It is possible in SLOP to implement this using the appropriate base class. In this example we will be implementing EvtGen's HQET2 for $B^0 \to D^{*}$ (for illustrative purposes as there is already CLN in SLOP). 

To make a FF scheme that works with SLOP, it must be a class that inherits from `FormFactor` class and implements the `get_ff` method. For $B^0 \to D^{*}$ we could inherit from `FormFactorBToV` but for full generality in this example we will be inheriting from the base `FormFactor` class.

```{code-block} python

from math import sqrt
from slophep.Predictions.FormFactorBase import FormFactor

class MyHQET2(FormFactor):
    def __init__(self, par: dict = None, scale: float = None):
        super().__init__(par, scale)

        self._name = "BToDst_HQET2"
        self._ffpar = {
            "RhoSq" : 1.122,
            "h_A1"  : 0.908,
            "R1"    : 1.270,
            "R2"    : 0.852,
            "R0"    : 1.15
        }
        self._params = ["RhoSq", "h_A1", "R1", "R2", "R0"]

        internalparams = {
            "Mb" : self.par["m_B0"],
            "Mc" : self.par["m_D*+"]
        }
        self.internalparams.update(internalparams)

    def get_ff(self, q2: float) -> dict:
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        w = (mB**2 + mV**2 - q2) / (2*mB*mV)
        z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2))
        RV = 2*sqrt(mB*mV)/(mB+mV)

        hA1_1 = self.ffpar['h_A1']
        R1_1 = self.ffpar['R1']
        R2_1 = self.ffpar['R2']
        R0_1 = self.ffpar['R0']
        rho2 = self.ffpar['RhoSq']

        hA1 = hA1_1 * (1 - 8*rho2*z + (53*rho2-15)*z**2 - (231*rho2-91)*z**3)
        R1 = R1_1 - 0.12*(w-1) + 0.05*(w-1)**2
        R2 = R2_1 + 0.11*(w-1) - 0.06*(w-1)**2
        R0 = R0_1 - 0.11*(w-1) + 0.01*(w-1)**2

        ff = {}
        ff['A1'] = hA1 * RV * (w+1)/2.
        ff['A0'] = R0/RV * hA1
        A2 = R2/RV * hA1
        # conversion from A_1, A_2 to A_12
        ff['A12'] = ((ff['A1']*(mB + mV)**2* (mB**2 - mV**2 - q2)
                - A2*(mB**4 + (mV**2 - q2)**2 - 2*mB**2*(mV**2 + q2)))
                / (16.*mB*mV**2*(mB + mV)))
        ff['V'] = R1/RV * hA1
        # SM only non tensor FFs
        ff["T1"] = 0.0
        ff["T2"] = 0.0
        ff["T23"] = 0.0

        return ff
```

In the example above, we set on initialisation the relevant FF parameters ($h_{A_1}$, $\rho^2$, $R_0(1)$, $R_1(1)$, and $R_2(1)$) such that they may be accessed through the `ffpar` property, as well as other relevant values to be accessed through the `internalparams` property.

After that it is a matter of implementing `MyHQET.get_ff` which returns a dictionary of the form-factors at a particular $q^2$. It is important here that the basis of the FFs must be $A_0$, $A_1$, $A_{12}$, $V$, $T_1$, $T_2$, $T_{23}$ as that is the basis SLOP uses to compute the amplitudes/observables.

Having implemented this class, we can now use it in predictions:
```{code-block} python

from slophep.Predictions.Observables import BdToDstEllNuPrediction

pred = BdToDstEllNuPrediction("mu", "mu", MyHQET2)
print(pred.J(q2=5.0))
```
which outputs:
```
{'1s': 0.37389761978504465, '1c': 0.5022136300317644, '2s': 0.12426178670587891, '2c': -0.49849696460619636, '6s': -0.3596762090456695, '6c': 0.0036411611548686847, 3: -0.1719490703269733, 4: 0.32373163006084216, 5: -0.27810535959217536, 7: -0.0, 8: 0.0, 9: -0.0}
```