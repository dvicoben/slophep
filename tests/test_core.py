import pytest

import slophep.Core.Parameter as spar
import slophep.Core.physdata as sphys

def test_parameter():
    p = spar.Parameter("testparam", 2.)
    assert p.name == "testparam"
    assert p.get_val() == 2.

    p.set_val(-1.)
    assert p.get_val() == -1.


def test_param_manager():
    pnames = ["m_B0", "m_D1*0"]
    pm = spar.ParameterManager()
    for pname in pnames:
        # Check that loaded defaults:
        assert pname in pm
        # Check loaded appropiate value
        assert pm[pname] == sphys.DEFAULT_PARAMS[pname]

    p = spar.Parameter("testparam", 2.)
    pm.add_param(p)
    # Test correctly retrieve param
    assert pm.get_param("testparam") is p

    # Test Aliases
    pm.add_alias("testparam", "testalias")
    assert "testalias" in pm
    # Test correctly retrieve param from alias
    assert pm.get_param("testalias") is p
    # Test remove alias
    pm.remove_alias("testalias")
    assert "testalias" not in pm

    # Test param registering
    pm.register_param("testparam2", 4.)
    assert "testparam2" in pm
    assert pm.get_val("testparam2") == 4.

    # Test not overwriting by default
    p2 = spar.Parameter("testparam", 3.)
    pm.add_param(p2)
    assert pm.get_val("testparam") != 3.
    # Test overwriting
    pm.add_param(p2, override=True)
    assert pm.get_val("testparam") == 3.
    assert pm.get_param("testparam") is p2

    # Test value setting
    pm.set_val("testparam", -1.)
    assert pm.get_val("testparam") == -1.
    # Test multiple value setting
    pm.set_vals({
        "testparam"  : -3.,
        "testparam2" : -4.
    })
    assert pm.get_val("testparam") == -3.
    assert pm.get_val("testparam2") == -4.

    # Fail to set value for inexistent param
    with pytest.raises(KeyError):
        pm.set_val("dummyparam", 2.)

    # Test WCs
    wcstr = "CVR_bcmunumu"
    wcoeffs = {
        wcstr : 0.1+0.5j,
    }
    pm.set_wc(wcoeffs)
    # Correctly register and set new WCs
    assert pm.get_val("WCRe:CVR_bcmunumu") == 0.1
    assert pm.get_val("WCIm:CVR_bcmunumu") == 0.5
    # Retrieve correct WCs
    assert pm._get_wc_values()[wcstr] == wcoeffs[wcstr]
    # Properly change existing WCs
    pm.set_wc({wcstr : 0.0})
    assert pm.get_val("WCRe:CVR_bcmunumu") == 0.0
    assert pm.get_val("WCIm:CVR_bcmunumu") == 0.0
    assert pm._get_wc_values()[wcstr] == 0.0


def test_param_manager_merge_into():
    pm1 = spar.ParameterManager()
    pm2 = spar.ParameterManager()
    pm1.register_param("test", 1.)
    pm2.register_param("test", 2.)

    pm_merged = pm1.merge_into(pm2)
    # Test we obtain the value from pm2
    assert pm_merged.get_val("test") == 2.


def test_param_manager_merge_onto():
    pm1 = spar.ParameterManager()
    pm2 = spar.ParameterManager()
    pm1.register_param("test", 1.)
    pm2.register_param("test", 2.)

    pm_merged = pm1.merge_onto(pm2)
    # Test we obtain the value from pm1
    assert pm_merged.get_val("test") == 1.


def test_param_user():
    pu = spar.ParameterUser()

    # Check getting method works
    assert pu.get_param("m_B0") == sphys.DEFAULT_PARAMS["m_B0"]

    # Change param manager
    pm = spar.ParameterManager()
    pu.set_parammanager(pm)
    assert pu.pm is pm

    # Test parameter setting
    pu.set_param("m_B0", 2.)
    assert pu.get_param("m_B0") == 2.