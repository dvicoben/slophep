"""
Some tests for FormFactor.
Note this checks functionality, NOT that the correct numbers are computed!
Likely to eventually add separate tests for different FFs a la https://github.com/dvicoben/slophep-checks
"""
import pytest
from slophep.Core.Parameter import ParameterManager
from slophep.Core.user_registry import FFregistry
from slophep.FormFactors.FormFactorBase import FormFactor
from slophep.Tools.errfluct_tools import FluctType

@pytest.fixture
def ff_instances() -> dict[str, FormFactor]:
    return FFregistry.get_instance_registry()

def test_ff_init(ff_instances: dict[str, FormFactor]):
    # Test that user parameters are produced
    for iffobj in ff_instances.values():
        idefaults = iffobj.user_params_defaults()
        assert type(idefaults) is dict
        assert len(idefaults) > 0


def test_ff_pm(ff_instances: dict[str, FormFactor]):
    # Test that values registered into pm if it is changed
    for iffobj in ff_instances.values():
        idefaults = iffobj.user_params_defaults()
        pm = ParameterManager()
        iffobj.set_parammanager(pm)
        assert iffobj.pm is pm
        # Now test it contains all the userparams it should
        for ipar in idefaults:
            assert ipar in iffobj.pm


def test_ff_setter(ff_instances: dict[str, FormFactor]):
    for iffobj in ff_instances.values():
        # Test setting from fully qualified name
        idefaults = iffobj.user_params_defaults()
        iparname = list(idefaults.keys())[0]
        iffobj.set_param(iparname, -99)
        assert iffobj.get_param(iparname) == -99

        # Test setting from userparam
        idefraw = iffobj.define_userparams()
        jparname = list(idefraw.keys())[0]
        iffobj.set_userparam(jparname, -90)
        assert iffobj.get_userparam(jparname) == -90

        # Test setting from set_ff
        iffobj.set_ff({jparname : -1})
        assert iffobj.get_userparam(jparname) == -1

        # And that things match fully qualified param
        assert iffobj.get_param(f"{iffobj.name}:{jparname}") == -1
        assert iffobj.get_param(f"{iffobj.name}:{jparname}") is iffobj.get_userparam(jparname)


def test_ff_compute(ff_instances: dict[str, FormFactor]):
    for iffobj in ff_instances.values():
        iffvals = iffobj.get_ff(5.0)
        # Test that proper output type
        assert type(iffvals) is dict
        assert len(iffvals) > 0
        # Check for a non-zero FF
        assert any([k != 0.0 for k in iffvals.values()])


# Might move this test to errorsampler test file
def test_ff_has_flucttype(ff_instances: dict[str, FormFactor]):
    for iffobj in ff_instances.values():
        # Check that fluctsettings decorator properly used
        assert hasattr(iffobj.get_ff, "_fluct_type")
        assert getattr(iffobj.get_ff, "_fluct_type") == FluctType.DICTNUMERIC