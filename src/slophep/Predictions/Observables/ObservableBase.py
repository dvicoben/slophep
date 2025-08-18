import flavio
from slophep.Predictions.FormFactorBase import FormFactor

class ObservableBase:
    def __init__(self, 
                 FF: FormFactor,
                 ffargs: list = [],
                 par: dict = None,
                 scale: float = 4.8,
                 ):
        """Base class for predictions to handle WC and FF setters"""
        self._par: dict = flavio.default_parameters.get_central_all()
        if type(par) == dict:
            self._par = par
        self._scale: float = scale

        self._FF: FormFactor = FF(par, scale, *ffargs)
        self._wc_obj: flavio.WilsonCoefficients = flavio.WilsonCoefficients()

    @property
    def scale(self) -> float: 
        """Renorm scale"""
        return self._scale
    @property
    def par(self) -> dict: 
        """Dictionary of parameters, defaults to `flavio.default_parameters.get_central_all()`"""
        return self._par
    @property
    def wc_obj(self) -> flavio.WilsonCoefficients: 
        """Wilson coefficient object"""
        return self._wc_obj
    @property
    def FF(self) -> FormFactor: 
        """Form factors"""
        return self._FF

    def set_wc(self, wc_dict: dict, eft: str = 'WET', basis: str = 'flavio'):
        """Set the wilson coefficients

        Parameters
        ----------
        wc_dict : dict
            Dictionary of wilson coefficients
        eft : str, optional
            EFT, by default 'WET'
        basis : str, optional
            WC basis (see https://wcxf.github.io/bases.html), by default 'flavio'
        """
        self._wc_obj.set_initial(wc_dict, self.scale, eft, basis)

    def set_ff(self, ffparams: dict):
        """Set form factor parameters. Can use None to leave a parameter unchanged.

        Parameters
        ----------
        ffparams : dict
            Dictionary of form factor parameters, names should match FF.ffpars dictionary in the particular 
            scheme being used
        """
        self._FF.set_ff(**ffparams)
    
    def _fullsetter(self, params: dict, constants: dict = {}):
        """Set WCs and FF parameters, for usage with Fluctuate

        FF parameters must follow naming in self.FF.params

        WCs must be in the flavio basis, prepended with 'WCRe\_' or 'WCIm\_'
        for the respective component

        Parameters
        ----------
        params : dict
            Dictionary of parameters to set
        constants: dict
            Dictinoary of parameters that are set constant - for specific use-cases with fluctuate
            e.g. in case want to set a particular WC to some non-zero value for all fluctations. In
            principle FF that are not in params should be kept the same so shouldn't need to
            pass them here.
        """
        # Note no key checking performed! User must make sure these are set-up correctly!
        par = {**params, **constants}
        wc = {}
        ff = {}
        for ikey, ival in par.items():
            # Handle Wilson Coefficients
            if "WCRe_" in ikey or "WCIm_" in ikey:
                name_wc = ikey[5:]
                iwc = ival if "WCRe_" in ikey else 1.0j*ival
                if name_wc not in wc:
                    wc[name_wc] = 0.0
                wc[name_wc] += iwc
            # Everything else is considered a FF
            else:
                ff[ikey] = ival
        
        if len(wc) > 0:
            self.set_wc(wc)
        if len(ff) > 0:
            self.set_ff(ff)


    def set_ff_fromlist(self, ffparams: list):
        """Set form factor parameters in form of list. Must include all parameters in self.FF.params, in order.
        Can use None to leave a parameter unchanged.

        Parameters
        ----------
        ffparams : list
            All FF parameters, in appropiate order
        """
        self._FF.set_ff_fromlist(ffparams)
