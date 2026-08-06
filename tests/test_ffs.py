"""
Some tests for FormFactor.
Note this checks functionality, NOT that the correct numbers are computed!
"""
import pytest
from slophep.Core.user_registry import FFregistry

@pytest.fixture
def ff_instances() -> dict:
    return FFregistry.get_instance_registry()

def test_ff_init(ff_instances):
    # Test that user parameters are produced correctly:
    for iffobj in ff_instances.values():
        idefaults = iffobj.user_params_defaults()
        assert type(idefaults) is dict
        assert len(idefaults) > 0

def test_ff_compute(ff_instances):
    for iffobj in ff_instances.values():
        iffvals = iffobj.get_ff(5.0)
        # Test that proper output type:
        assert type(iffvals) is dict
        assert len(iffvals) > 0
        # Check for a non-zero FF
        assert any([k != 0.0 for k in iffvals.values()])