"""
Some tests for Observables. 
A bit bare-bones at the moment, as each mode may require bespoke test.
Note this checks functionality, NOT that the correct numbers are computed!
"""
import slophep.Observables as sobs
import slophep.FormFactors as sff
from slophep.Core.Parameter import ParameterManager
from slophep.Observables.ObservableImpl.ObservableBase import Observable
from slophep.Tools.errfluct_tools import FluctType

# A bit sloppy, probably best to also make some kind of registry...
def obs_instances_withang(l) -> dict[str, Observable]:
    obs = [
        sobs.BdToDstEllNuPrediction(l, l, sff.BdToDstFF.BLPR()),
        sobs.BsToDsstEllNuPrediction(l, l, sff.BsToDsstFF.BLPR()),
        sobs.BdToDEllNuPrediction(l, l, sff.BdToDFF.BLPR()),
        sobs.BuToDEllNuPrediction(l, l, sff.BuToDFF.BLPR()),
        sobs.BdToPiEllNuPrediction(l, l, sff.BdToPiFF.BSZ()),
        sobs.LbToLcEllNuPrediction(l, l, sff.LbToLcFF.DKMR()),
    ]
    return {iobs.name : iobs for iobs in obs}

def obs_instances_noang(l) -> dict[str, Observable]:
    obs = [
        sobs.BcToJpsiEllNuPrediction(l, l, sff.BcToJpsiFF.HPQCD2020()),
        sobs.BsToKEllNuPrediction(l, l, sff.BsToKFF.BCL()),
        sobs.BdToD0stEllNuPrediction(l, l, sff.BdToD0stFF.ISGW2()), 
        sobs.BdToD1EllNuPrediction(l, l, sff.BdToD1FF.ISGW2()), 
        sobs.BdToD1stEllNuPrediction(l, l, sff.BdToD1stFF.ISGW2()), 
        sobs.BdToD2stEllNuPrediction(l, l, sff.BdToD2stFF.ISGW2()),
        sobs.BuToD0stEllNuPrediction(l, l, sff.BuToD0stFF.ISGW2()), 
        sobs.BuToD1EllNuPrediction(l, l, sff.BuToD1FF.ISGW2()), 
        sobs.BuToD1stEllNuPrediction(l, l, sff.BuToD1stFF.ISGW2()), 
        sobs.BuToD2stEllNuPrediction(l, l, sff.BuToD2stFF.ISGW2()),
        sobs.BsToDs0stEllNuPrediction(l, l, sff.BsToDs0stFF.ISGW2()), 
        sobs.BsToDs1EllNuPrediction(l, l, sff.BsToDs1FF.ISGW2()), 
        sobs.BsToDs1stEllNuPrediction(l, l, sff.BsToDs1stFF.ISGW2()), 
        sobs.BsToDs2stEllNuPrediction(l, l, sff.BsToDs2stFF.ISGW2()),
    ]
    return {iobs.name : iobs for iobs in obs}

def obs_instances_all(l) -> dict[str, Observable]:
    obs = {}
    obs.update(obs_instances_noang(l))
    obs.update(obs_instances_withang(l))
    return obs



def test_obs_pm():
    # Tests the parameter manager is shared with internal FFs
    obs = obs_instances_all("mu")
    for iname, iobs in obs.items():
        # Test pm is the same
        assert iobs.pm is iobs.FF.pm
        # Test setting a new one
        pm = ParameterManager()
        iobs.set_parammanager(pm)
        assert iobs.pm is pm
        assert iobs.pm is iobs.FF.pm


def test_obs_dG():
    obs = obs_instances_all("mu")
    query_points = [2., 6., 9.]
    for iname, iobs in obs.items():
        # Test we get some prediction, several points queried
        qpoints = [iobs.dGdq2(iqsq) for iqsq in query_points]
        # Test rate positive at all queried points
        assert all([k >= 0.0 for k in qpoints])
        # Test we get non-zero somewhere
        assert any([k > 0.0 for k in qpoints])


def test_obs_ang():
    obs = obs_instances_withang("mu")
    for iname, iobs in obs.items():
        # Test we get some prediction, several points queried
        ipred = iobs.J(8.0)
        assert type(ipred) is dict
        assert len(ipred) > 0
        # Could check for a non-zero but may break depending on q2 range...


def test_dG_has_flucttype():
    obs = obs_instances_all("mu")
    for iobj in obs.values():
        # Check that fluctsettings decorator properly used
        assert hasattr(iobj.dGdq2, "_fluct_type")
        assert getattr(iobj.dGdq2, "_fluct_type") == FluctType.NUMERIC


def test_J_has_flucttype():
    obs = obs_instances_withang("mu")
    for iobj in obs.values():
        # Check that fluctsettings decorator properly used
        assert hasattr(iobj.J, "_fluct_type")
        assert getattr(iobj.J, "_fluct_type") == FluctType.DICTNUMERIC