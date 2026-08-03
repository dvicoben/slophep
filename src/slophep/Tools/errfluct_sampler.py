from __future__ import annotations
import numpy as np
from enum import StrEnum
import json

import slophep.Tools.errfluct_tools as errflt
from slophep.Core.Parameter import ParameterUser

import logging
logger = logging.getLogger(__name__)



class ErrorType(StrEnum):
    SYMMETRIC        = "Symmetric" ,
    GAUSSIAN         = "Symmetric" ,
    ASYMMETRIC       = "Asymmetric",
    BIFURCATED_GAUSS = "Asymmetric",
    CUSTOM           = "Custom"    



class ErrorSampler:
    def __init__(self, 
                 params     : list[str], 
                 mean       : dict[str, float],
                 fluctuator : errflt.Fluctuator,
                 constants  : dict[str, float] = {}):
        self._fluctuations  : list              = []
        self._params        : list[str]         = params
        self._mean          : dict[str, float]  = mean
        self._constants     : dict[str, float]  = constants
        self._fluctuator    : errflt.Fluctuator = fluctuator

    @property
    def fluctuations(self) -> list:
        return self._fluctuations
    @property
    def params(self) -> list[str]:
        return self._params
    @property
    def fluctuator(self) -> errflt.Fluctuator:
        return self._fluctuator
    @property
    def mean(self) -> dict[str, float]:
        return self._mean
    @property
    def constants(self) -> dict[str, float]:
        return self._constants

    def clear_fluctuations(self):
        """Clear produced fluctuations"""
        self._fluctuations = []

    def fluctuate(self, 
                  N: int, 
                  rng: np.random.Generator = np.random.default_rng()) -> None:
        """Produce N gaussian/bifurcated gaussian fluctuations

        Parameters
        ----------
        N : int
            Number of fluctuations
        rng: np.random.Generator, optional
            numpy generator
        """
        flucts = self.fluctuator.generate(N, rng)
        if len(flucts) < 1:
            logger.info("No fluctuations generated, no changes")
            return

        self._fluctuations = np.array(flucts)

    def compute_fluctations(self, 
                            param_user : ParameterUser, 
                            attr       : str, 
                            attr_args  : list = []) -> list:
        logger.info(f"Computing {attr} with {attr_args} for {len(self.fluctuations)} fluctuations")
        res = []
        init_vals = {param_user.get_param(kpar) for kpar in self.mean}
        init_vals.update({param_user.get_param(kpar) for kpar in self.constants})

        for ifluct in self.fluctuations:
            ipar = {self.params[k] : ifluct[k] for k in range(len(self.params))}
            ipar.update(self.constants)
            # print(ipar)
            param_user.pm.set_vals(ipar)
            feval = getattr(param_user, attr)
            x = feval(*attr_args) if len(attr_args) > 0 else feval()
            res.append(x)
        
        param_user.pm.set_vals(init_vals)

        return res

    def central_val(self, 
                    param_user : ParameterUser, 
                    attr       : str, 
                    attr_args  : list = []):
        init_vals = {param_user.get_param(kpar) for kpar in self.mean}
        init_vals.update({param_user.get_param(kpar) for kpar in self.constants})
        vals = {**self.mean, **self.constants}
        param_user.pm.set_vals(vals)
        feval = getattr(param_user, attr)
        x = feval(*attr_args) if len(attr_args) > 0 else feval()
        param_user.pm.set_vals(init_vals)
        return x

    def get_error(self, 
                  param_user : ParameterUser, 
                  attr       : str, 
                  attr_args  : list = [], 
                  cl         : float = 0.683):
        """Compute error bands using sampling for a particular user method

        Parameters
        ----------
        param_user : ParameterUser
            The ParameterUser whose values will be fluctuated
        attr : str
            The attribute/method in self.obs to compute
        attr_args : list
            Arguments to be passed to that attribute, by default an empty
            list which means no arguments to pass
        cl : float, optional
            CL to asses error for, by default 0.683
        """ 
        res = self.compute_fluctations(param_user, attr, attr_args)
        alpha = 1.0-cl
        feval = getattr(param_user, attr)

        if not hasattr(feval, "_fluct_type"):
            logger.warning(f"ParameterUser {param_user}.{attr} has no defined _fluct_type, returning all fluctuations")
            return res

        fluct_type = getattr(feval, "_fluct_type")
        return errflt._FLUCTTYPE_ERRORCALC[fluct_type](res, alpha)


    # Factory methods - should be preferred way of setting-up sampler
    @classmethod
    def create_custom(cls, 
            params       : list[str], 
            mean         : dict[str, float],
            fluctuations : list = [],
            constants    : dict[str, float] = {}) -> ErrorSampler:
        sampler = cls(params, mean, errflt.Fluctuator(), constants)
        sampler._fluctuations = fluctuations
        return sampler
    
    @classmethod
    def create_symmetric_from_cov(cls,
            params    : list[str],
            mean      : dict[str, float],
            cov       : np.ndarray,
            constants : dict[str, float] = {}) -> ErrorSampler:
        mvec = [mean[ipar] for ipar in params]
        fluctuator = errflt.FlucutatorGaussian(mvec, cov)
        sampler = cls(params, mean, fluctuator, constants)
        return sampler

    @classmethod
    def create_symmetric_from_corr(cls,
            params    : list[str],
            mean      : dict[str, float],
            err       : dict[str, float],
            corr      : np.ndarray,
            constants : dict[str, float] = {}) -> ErrorSampler:
        mvec = [mean[ipar] for ipar in params]
        errvec = [err[ipar] for ipar in params]
        cov = errflt.cov_from_corr(corr, np.array(errvec))
        fluctuator = errflt.FlucutatorGaussian(mvec, cov)
        sampler = cls(params, mean, fluctuator, constants)
        return sampler
                
    @classmethod
    def create_asymmetric(cls,
            params    : list[str],
            mean      : dict[str, float],
            errlo     : dict[str, float],
            errhi     : dict[str, float],
            corr      : np.ndarray,
            constants : dict[str, float] = {}) -> ErrorSampler:
        mvec = [mean[ipar] for ipar in params]
        errlovec = [errlo[ipar] for ipar in params]
        errhivec = [errhi[ipar] for ipar in params]
        fluctuator = errflt.FlucutatorBifurGaussian(mvec, errlovec, errhivec, np.array(corr))
        sampler = cls(params, mean, fluctuator, constants)
        return sampler

    @classmethod
    def create_from_configfile(cls, inpath: str) -> ErrorSampler:
        with open(inpath) as f:
            config: dict = json.load(f)

        params = config["Params"]
        mean = config["Mean"]        
        constants = config.get("Constants", {})
        errtype_str = config.get("Type")
        if errtype_str is None:
            logger.info("Unspecified error Type, assuming Symmetric")
            errtype_str = ErrorType.SYMMETRIC
        errtype = ErrorType(errtype_str)

        sampler = None
        match errtype:
            case ErrorType.SYMMETRIC:
                if config.get("Cov") is not None:
                    sampler = cls.create_symmetric_from_cov(params, mean, config["Cov"], constants)
                else:
                    sampler = cls.create_symmetric_from_corr(params, mean, config["Error"], config["Cov"], constants)
            case ErrorType.ASYMMETRIC:
                sampler = cls.create_asymmetric(params, mean, config["Error_lo"], config["Error_hi"], config["Corr"])
            case ErrorType.CUSTOM:
                sampler = cls.create_custom(params, mean, config["Fluctuations"], constants)
        
        return sampler

    @classmethod
    def create_blank_config(cls, outpath: str, errtype: ErrorType = ErrorType.SYMMETRIC) -> dict:
        config = {
            "Comments"  : [],
            "Type"      : str(errtype),
            "Params"    : [],
            "Mean"      : {},
            "Constants" : {},
        }
        match errtype:
            case ErrorType.SYMMETRIC:
                config.update({
                    "Cov" : []
                })
            case ErrorType.ASYMMETRIC:
                config.update({
                    "Error_lo" : {},
                    "Error_hi" : {},
                    "Corr" : []
                })
            case ErrorType.CUSTOM:
                config.update({
                    "Fluctuations" : []
                })

        with open(outpath, "w") as outfile:
            json.dump(config, outfile, indent=2)
        return config
        