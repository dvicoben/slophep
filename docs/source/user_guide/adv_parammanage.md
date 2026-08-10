# Parameters in SLOP

For some of the more advanced use, like implementing a custom form-factor scheme, it is useful to be familiar with how parameters are managed in SLOP. A short summary is provided here.

## Parameter Managers

Parameters are handled by a `ParameterManager`. Parameters are added/registered to the manager, and then used by a `ParameterUser` to compute quantities of interest. This is done to ensure that different `ParameterUser`s can use the same parameters simply by drawing them from the same `ParameterManager`, rather than having to ensure that different internal attributes of each user are all aligned.

Parameters in the `ParameterManager` must have distinct names. By default, a `ParameterManager` will load in a series of physics parameters (e.g. masses, lifetimes), which you can print out:
```{code-block} python

from slophep.Core.Parameter import Parameter, ParameterManager

pm = ParameterManager()
pm.print_values()
```
To ensure parameters are not accidentally overriden, by default you cannot add/register parameters already in the manager. For example, try running the following:
```{code-block} python

print("Value before attempting to add parameter:", pm.get_val("m_B0"))
my_B0_M = Parameter("m_B0", 2.)
pm.add_param(my_B0_M)
print("Value after attempting to add parameter:", pm.get_val("m_B0"))
pm.register_param("m_B0", 2.)
print("Value after attempting to register parameter:", pm.get_val("m_B0"))
```
In none of the two attempts did the value change. We can change the value using `pm.set_val("m_B0", 2.)`, but not replace it by accident. It can be replaced intentionally using `pm.add_param(my_B0_M, override=True)` if it really is desired.


## Parameter Users

A `ParameterUser` is a base class that takes values from an internal `ParameterManager` and uses them for the computation of some quantity of interest. In SLOP, observables/predictions and form-factors are `ParameterUser`s. A `ParameterUser` should have a unique name that is class specific (not instance specific), and is stored in the class attribute `_name`. Because there are unique parameters for separate parameter users, classes inheriting from `ParameterUser` implement the `define_userparams()` method to indicate parameters that should be registered to the manager (if not already present) on initialization. 

As an example let's look at a form-factor class:
```{code-block} python

from slophep.FormFactors import BdToDstFF

cln = BdToDstFF.CLN()
# Now lets look at what's in the cln ParameterManager
cln.pm.print_values()
```
You will notice now there are a series of parameters at the end of the list:
```
  FFBdToDst@CLN:RhoSq                  : 1.207
  FFBdToDst@CLN:h_A1                   : 0.908
  FFBdToDst@CLN:R1                     : 1.401
  FFBdToDst@CLN:R2                     : 0.854
  FFBdToDst@CLN:R0                     : 1.15
```
which are prefixed with `FFBdToDst@CLN`. Let's look now at the output of `cln.define_userparams()`, which is:
```
{'RhoSq': 1.207, 'h_A1': 0.908, 'R1': 1.401, 'R2': 0.854, 'R0': 1.15}
```
On initialization, a `ParameterUser` takes the output of `define_userparams()`, and registers everything in it, prefixed with the name of the `ParameterUser`. The convention is that the naming of parameters is prefixed with the name of the `ParameterUser`, i.e. of the form `{NAME_OF_ParameterUser}:{NAME_OF_PARAMETER}` to avoid naming collisions with all loaded default physics parameters (e.g. masses, lifetimes), which do not have a prefix. 

A `ParameterUser` can access parameters in its internal `ParameterManager` through `get_param` or `get_userparam`:
```{code-block} python

print(cln.get_param("m_B0"))
print(cln.get_param("FFBdToDst@CLN:RhoSq"))
print(cln.get_userparam("RhoSq"))
```
Using `get_param` requires the full, qualified name. The method `get_userparam` is simply a convenient short-cut to query user-specific parameters without having to input the `ParameterUser`'s name (under the hood it is simply `get_userparam(param)` is simply `get_param(f"{self.name}:{param}")`). One may change parameter values similarly through `set_param` and `set_userparam`.

```{note}

You can access a `ParameterUser`'s internal `ParameterManager` through it's `pm` property. However, one should not, in general, interface with the user's parameters like this, but through the `ParameterUser` methods `get_param`, `set_param`, `get_userparam` and `set_userparam`.
```

You can change the manager of a `ParameterUser` through `set_parammanager`:
```{code-block} python

cln.set_parammanger(pm)
```
which also takes care to register any absent parameters required by the user into the new manager.


## About the `ErrorSampler`

The `ErrorSampler` works by interacting directly with the `ParameterManager` of the `ParameterUser` object whose methods we want to sample. Because of this, names of fluctuating parameters supplied to the `ErrorSampler`, and in the `json` config-files for it, need to be the fully qualified names.