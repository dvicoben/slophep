import numpy as np
import json

class FluctuatePrediction:
    def __init__(self, obs):
        self._params = []
        self._cov = None
        self._mean = {}
        self._fluctuations = []
        self._obs_obj = obs
        self._constants = {}

    @property
    def params(self) -> list[str]:
        """List of fluctuated params names, in order as they appear in cov matrix"""
        return self._params
    @property
    def cov(self) -> np.ndarray:
        """Covariance matrix of parameters"""
        return self._cov
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
        return self._obs_obj

    def set_params(self, param_names: list[str], mean: dict, cov: list | np.ndarray, constants: dict = {}):
        """Sets parameters necessary for fluctuation

        Parameters
        ----------
        param_names : list[str]
            "List of fluctuated params names, in order as they appear in cov matrix
        mean : dict
            Nominal/mean values for the parameters
        cov : list | np.ndarray
            Covariance matrix
        """
        self._mean = mean
        self._params = param_names
        self._cov = np.array(cov)
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
        self.set_params(config["Params"], config["Mean"], config["Cov"], constants)

    def clear_fluctuations(self):
        """Clear produced fluctuations"""
        self._fluctuations = []

    def fluctuate(self, N: int, seed: int = None):
        rng = np.random.default_rng(seed)
        mean = [self.mean[k] for k in self.params]
        fluct = rng.multivariate_normal(mean, self.cov, int(N))
        self._fluctuations = fluct
    
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
            CL to asses error for, by default 0.68
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
            self.obs._fullsetter(ipar)
            feval = getattr(self.obs, attr)
            x = feval(*attr_args) if len(attr_args) > 0 else feval()
            res.append(x)
        
        self.obs._fullsetter(self.mean)

        if return_all:
            return res

        # Beware this breaking if the entires in the dictionary are not numeric
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
            return (res[int(len(res)*0.5*cl)], res[int(len(res)*(1-0.5*cl))])
        except:
            print("Fluctuations requested don't have support for anything other than return_all, returning all")
            return res
        
        
        

