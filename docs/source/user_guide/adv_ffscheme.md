# Custom FF Schemes

It may be that there is a particular form-factor scheme you want to implement. It is possible in SLOP to implement this using the appropriate base class. In this example we will be implementing EvtGen's HQET2 for $B^0 \to D^{*}$ (for illustrative purposes as there is already CLN in SLOP). 

## Making the FF Scheme

To make a FF scheme that works with SLOP, it must be a class that inherits from `FormFactor` class and implements the `get_ff` method. For $B^0 \to D^{*}$ we could inherit from `FormFactorBToV` but for full generality in this example we will be inheriting from the base `FormFactor` class.

```{code-block} python

from math import sqrt
from slophep.FormFactors import FormFactor

class MyHQET2(FormFactor):
    _name = "FFMyHQET"
    def define_userparams(self) -> dict[str, Any]:
        ffpar = {
            "RhoSq" : 1.122,
            "h_A1"  : 0.908,
            "R1"    : 1.270,
            "R2"    : 0.852,
            "R0"    : 1.15
        }
        return ffpar

    def get_ff(self, q2: float) -> dict[str, float]:
        mB = self.get_param("m_B0")
        mV = self.get_param("m_D*+")
        w = (mB**2 + mV**2 - q2) / (2*mB*mV)
        z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2))
        RV = 2*sqrt(mB*mV)/(mB+mV)

        hA1_1 = self.get_userparam('h_A1')
        R1_1  = self.get_userparam('R1')
        R2_1  = self.get_userparam('R2')
        R0_1  = self.get_userparam('R0')
        rho2  = self.get_userparam('RhoSq')

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
We give our form-factor implementation, `MyHQET` a unique name `_name = "FFMyHQET"`. Then, we use the method `define_userparams` to register the FF parameters ($h_{A_1}$, $\rho^2$, $R_0(1)$, $R_1(1)$, and $R_2(1)$). They will be registered in an internal `ParameterManager`, with fully qualified names the form `{_name}:{parameter}`, e.g. in this case `MyHQET2:RhoSq`, `MyHQET2:h_A1`, etc. Parameters defined here should be accessed through the `get_userparam` to be used elsewhere (e.g. `self.get_userparam("RhoSq")`), or through `get_param` using the fully qualified name, (e.g. `self.get_param("MyHQET2:RhoSq")`). Other relevant values, such as meson masses, will likely be in the class's internal `ParameterManager` already (if not you can also define them in `define_userparams`), and you can access them as using `get_param` as in the example.

After that it is a matter of implementing `MyHQET.get_ff` which returns a dictionary of the form-factors at a particular $q^2$. It is important here that the basis of the FFs must be $A_0$, $A_1$, $A_{12}$, $V$, $T_1$, $T_2$, $T_{23}$ as that is the basis SLOP uses to compute the amplitudes/observables.

Having implemented this class, we can now use it in predictions:
```{code-block} python

from slophep.Observables import BdToDstEllNuPrediction

pred = BdToDstEllNuPrediction("mu", "mu", MyHQET2())
print(pred.J(q2=5.0))
```
which outputs:
```
{'1s': 0.37389761978504465, '1c': 0.5022136300317644, '2s': 0.12426178670587891, '2c': -0.49849696460619636, '6s': -0.3596762090456695, '6c': 0.0036411611548686847, 3: -0.1719490703269733, 4: 0.32373163006084216, 5: -0.27810535959217536, 7: -0.0, 8: 0.0, 9: -0.0}
```

## Working with the `ErrorSampler`

If we want to generate errorbands, we can make a few additions to interface more nicely with the `ErrorSampler`. We will want to import the `@fluctsettings` decorator to declare the kind of output we expect from `get_ff`

```{code-block} python

from math import sqrt
from slophep.FormFactors import FormFactor
from slophep.Tools.errfluct_tools import FluctType, fluctsettings # We import some necessary utilites

class MyHQET2(FormFactor):
    _name = "FFMyHQET"
    def define_userparams(self) -> dict[str, Any]:
        ffpar = {
            "RhoSq" : 1.122,
            "h_A1"  : 0.908,
            "R1"    : 1.270,
            "R2"    : 0.852,
            "R0"    : 1.15
        }
        return ffpar

    @fluctsettings(FluctType.DICTNUMERIC)  # Only need this minor addition here!
    def get_ff(self, q2: float) -> dict[str, float]:
        mB = self.get_param("m_B0")
        mV = self.get_param("m_D*+")
        w = (mB**2 + mV**2 - q2) / (2*mB*mV)
        z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2))
        RV = 2*sqrt(mB*mV)/(mB+mV)

        hA1_1 = self.get_userparam('h_A1')
        R1_1  = self.get_userparam('R1')
        R2_1  = self.get_userparam('R2')
        R0_1  = self.get_userparam('R0')
        rho2  = self.get_userparam('RhoSq')

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
Now the `MyHQET2.get_ff` method can interface with the `ErrorSampler` in a manner similar to the example in [Producing Errorbands](baseuse-producing-errorbands). This decorator ensures that the output of using the `ErrorSampler` on `get_ff` will be of the form `{"FFname" : (error_lo, error_hi)}`, i.e. an output of the type `dict[str, tuple[float, float]]`. Note that you can still use the `ErrorSampler` even without this edit, but it will simply return a list with the output of `get_ff` for every sample/fluctuation, instead of automatically processing it into the dictionary with lower and upper bounds.