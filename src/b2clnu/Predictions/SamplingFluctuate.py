import numpy as np
import json
from b2clnu.Predictions.SamplingTools import bifurcated_gaussian_sampler

class SamplingHelper:
    def __init__(self, obs):
        self._params = []
        self._errtype = "Symmetric"
        self._cov = None
        self._corr = None
        self._mean = {}
        self._err = {}
        self._fluctuations = []
        self._obs_obj = obs
        self._constants = {}

    @property
    def params(self) -> list[str]:
        """List of fluctuated params names, in order as they appear in cov matrix"""
        return self._params
    @property
    def errtype(self) -> str:
        """Error type, should be `Symmetric` or `Asymmetric`"""
        return self._errtype
    @property
    def err(self) -> dict[float]:
        """Error values for the parameters"""
        return self._err
    @property
    def cov(self) -> np.ndarray:
        """Covariance matrix of parameters"""
        return self._cov
    @property
    def corr(self) -> np.ndarray:
        """Correlation matrix of parameters"""
        return self._corr
    @property
    def mean(self) -> dict[float]:
        """Nominal/mean values for the parameters"""
        return self._mean
    @property
    def nominal(self) -> dict[float]:
        """Alias for mean"""
        return self._mean
    @property
    def constants(self) -> dict[float]:
        """Parameters that are kept constant in fluctuation and their values.
           If they are here, they SHOULD NOT be in self.params
        """
        return self._constants
    @property
    def fluctuations(self) -> list | np.ndarray:
        """List of fluctuations generated"""
        return self._fluctuations
    @property
    def obs(self):
        """Observable/FF object to fluctuate"""
        return self._obs_obj

    def set_params_symmetric(self, param_names: list[str], mean: dict, cov: list | np.ndarray, constants: dict = {}):
        """Sets parameters necessary for fluctuation, for symmetric gaussian errors

        Parameters
        ----------
        param_names : list[str]
            List of fluctuated params names, in order as they appear in cov matrix
        mean : dict
            Nominal/mean values for the parameters
        cov : list | np.ndarray
            Covariance matrix
        constants: dict, optional
            Dictinoary of parameters that are set constant - for specific use-cases with fluctuate
            e.g. in case want to set a particular WC to some non-zero value for all fluctations. In
            principle FF that are not in params should be kept the same so shouldn't need to
            pass them here.
        """
        self._errtype = "Symmetric"
        self._mean = mean
        self._params = param_names
        self._cov = np.array(cov)
        self._constants = constants

    def set_params_asymmetric(self, param_names: list[str], mean: dict, errlo: dict, errhi: dict, corr: list | np.ndarray, constants: dict = {}):
        """Sets parameters necessary for fluctuation, for asymmetric gaussian errors

        Parameters
        ----------
        param_names : list[str]
            "List of fluctuated params names, in order as they appear in cov matrix
        mean : dict
            Nominal/mean values for the parameters
        errlo : dict
            Size of lower error
        errhi : dict
            Size of upper error
        corr : list | np.ndarray
            Correlation matrix
        constants: dict, optional
            Dictinoary of parameters that are set constant - for specific use-cases with fluctuate
            e.g. in case want to set a particular WC to some non-zero value for all fluctations. In
            principle FF that are not in params should be kept the same so shouldn't need to
            pass them here.
        """
        self._errtype = "Asymmetric"
        self._params = param_names
        self._mean = mean
        self._err["lo"] = errlo
        self._err["hi"] = errhi
        self._corr = np.array(corr)
        self._constants = constants
    
    def set_params_from_configfile(self, inpath: str):
        """Set parameters for fluctuation from a json file

        Parameters
        ----------
        inpath : str
            Path to json file
        """
        with open(inpath) as f:
            config = json.load(f)
        constants = config["Constants"] if "Constants" in config else {}
        if config["Type"] == "Asymmetric":
            self.set_params_asymmetric(config["Params"], config["Mean"],
                                       config["Error_lo"], config["Error_hi"], config["Corr"],
                                       constants)
        elif config["Type"] == "Symmetric":
            self.set_params_symmetric(config["Params"], config["Mean"], config["Cov"], constants)
        else:
            print("Specified error type is not Asymmetric or Symmetric, assuming Symmetric")
            self.set_params_symmetric(config["Params"], config["Mean"], config["Cov"], constants)

    def clear_fluctuations(self):
        """Clear produced fluctuations"""
        self._fluctuations = []

    def fluctuate(self, N: int, seed: int = None):
        """Produce N gaussian/bifurcated gaussian fluctuations

        Parameters
        ----------
        N : int
            Number of fluctuations
        seed : int, optional
            numpy random seed, by default None
        """
        rng = np.random.default_rng(seed)
        mean = np.array([self.mean[k] for k in self.params])
        fluct = []
        if self.errtype == "Asymmetric":
            errlo = np.array([self.err["lo"][k] for k in self.params])
            errhi = np.array([self.err["hi"][k] for k in self.params])
            fluct = [bifurcated_gaussian_sampler(mean, errlo, errhi, self.corr, rng) for isample in range(N)]
        else:
            fluct = rng.multivariate_normal(mean, self.cov, int(N))
        
        self._fluctuations = np.array(fluct)
    
    def get_error(self, attr: str, attr_args: list = [], cl: float = 0.683, return_all: bool = False):
        """Gets an error using produced fluctations for a spcified CL

        Parameters
        ----------
        attr : str
            The attribute/method in self.obs to compute
        attr_args : list
            Arguments to be passed to that attribute, by default an empty
            list which means no arguments to pass
        cl : float, optional
            CL to asses error for, by default 0.683
        return_all : bool, optional
            Whether to return a list with all the fluctiations (only) 
            rather than the errors, by default False
        """
        print(f"Computing {attr} with {attr_args} for {len(self.fluctuations)} fluctuations")
        res = []
        alpha = 1.0-cl

        for ifluct in self.fluctuations:
            ipar = {self.params[k] : ifluct[k] for k in range(len(self.params))}
            # print(ipar)
            self.obs._fullsetter(ipar, self.constants)
            feval = getattr(self.obs, attr)
            x = feval(*attr_args) if len(attr_args) > 0 else feval()
            res.append(x)
        
        self.obs._fullsetter(self.mean, self.constants)

        if return_all:
            return res

        # Beware this breaking if the entries in the dictionary are not numeric
        if type(res[0]) == dict:
            merged_res = {k: np.sort([d.get(k) for d in res])
                          for k in res[0]} #set().union(*res)
            # Get upper and lower errors
            errs = {
                k : (ires[int(len(ires)*0.5*alpha)], ires[int(len(ires)*(1-0.5*alpha))]) 
                for k, ires in merged_res.items()
            }
            return errs
        
        # Use this to catch any kind of numeric types (numpy.float64, float, etc.)
        try:
            tmp = int(res[0])
            res = np.sort(res)
            return (res[int(len(res)*0.5*alpha)], res[int(len(res)*(1-0.5*alpha))])
        except:
            print("Fluctuations requested don't have support for anything other than return_all, returning all")
            return res
        
        
        

