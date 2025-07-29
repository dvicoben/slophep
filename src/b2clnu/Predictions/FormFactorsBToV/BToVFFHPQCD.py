from b2clnu.Predictions.FormFactorsBToV import FormFactorBToV
from math import sqrt


class HPQCD_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, V, par, scale)
        
        self._name = "HPQCD"
        # From supplementary material in https://arxiv.org/abs/2304.03137v2 
        # This goes up to order 10 in (w-1) for each FF !!
        # Will work on a neater way to store these values ...
        self._ffpar = { 
            "a^0_hA1"  : 0.9033127209857411,
            "a^1_hA1"  : -0.9508435406529809,
            "a^2_hA1"  : 0.5257520677535575,
            "a^3_hA1"  : -0.10086711625761521,
            "a^4_hA1"  : -0.18939192845292227,
            "a^5_hA1"  : 0.3400524945787585,
            "a^6_hA1"  : -0.38728694249597406,
            "a^7_hA1"  : 0.37132088878018327,
            "a^8_hA1"  : -0.3237217925179691,
            "a^9_hA1"  : 0.2653479192536277,
            "a^10_hA1" :  -0.20819351685286072,
            "a^0_hA2"  : -0.4642005293547826,
            "a^1_hA2"  : 0.4762715910680349,
            "a^2_hA2"  : -1.024703497027969,
            "a^3_hA2"  : 1.4647317764574475,
            "a^4_hA2"  : -1.655288306586653,
            "a^5_hA2"  : 1.6308234881825487,
            "a^6_hA2"  : -1.4706335052056245,
            "a^7_hA2"  : 1.2471707316030383,
            "a^8_hA2"  : -1.011253949245923,
            "a^9_hA2"  : 0.7925552531002823,
            "a^10_hA2" : -0.6049480665398556,
            "a^0_hA3"  : 1.3554401987030802,
            "a^1_hA3"  : -1.2307794387685247,
            "a^2_hA3"  : 0.4092337281757482,
            "a^3_hA3"  : 0.3128491447530095,
            "a^4_hA3"  : -0.7495739296382179,
            "a^5_hA3"  : 0.9267631598398023,
            "a^6_hA3"  : -0.9269882683956783,
            "a^7_hA3"  : 0.8292611271071676,
            "a^8_hA3"  : -0.6910768093514228,
            "a^9_hA3"  : 0.5479508577968928,
            "a^10_hA3" : -0.41856224634866757,
            "a^0_hV"   : 1.2231306553287695,
            "a^1_hV"   : -1.9958578437965673,
            "a^2_hV"   : 2.55917256361585,
            "a^3_hV"   : -2.671400245583292,
            "a^4_hV"   : 2.4971801805018696,
            "a^5_hV"   : -2.181267497854233,
            "a^6_hV"   : 1.8206590102598812,
            "a^7_hV"   : -1.471445852774884,
            "a^8_hV"   : 1.1612767702518605,
            "a^9_hV"   : -0.9001266450293959,
            "a^10_hV"  : 0.6880569014570732,
            "a^0_hT1"  : 0.9229006472316218,
            "a^1_hT1"  : -1.1546465073714087,
            "a^2_hT1"  : 0.9681015610543681,
            "a^3_hT1"  : 0.7695270465359492,
            "a^4_hT1"  : -0.525549938591732,
            "a^5_hT1"  : -0.3947275978604678,
            "a^6_hT1"  : -0.17881029770168894,
            "a^7_hT1"  : -0.06068235137612874,
            "a^8_hT1"  : -0.014820648861609841,
            "a^9_hT1"  : -0.0010156037902605876,
            "a^10_hT1" : 0.0017528617366402374,
            "a^0_hT2"  : -0.10972353521699603,
            "a^1_hT2"  : 0.029878371191936617,
            "a^2_hT2"  : -0.3914898707269266,
            "a^3_hT2"  : 7.435671398299036,
            "a^4_hT2"  : 4.1672064427294435,
            "a^5_hT2"  : 1.8337356857119789,
            "a^6_hT2"  : 0.7240526579747844,
            "a^7_hT2"  : 0.2710686655731731,
            "a^8_hT2"  : 0.09883480662543806,
            "a^9_hT2"  : 0.03557789930656164,
            "a^10_hT2" : 0.012731997593546481,
            "a^0_hT3"  : -0.20107689611999296,
            "a^1_hT3"  : 1.30505464841678,
            "a^2_hT3"  : -0.8206859343721933,
            "a^3_hT3"  : 1.6015122623914702,
            "a^4_hT3"  : 1.442334286115229,
            "a^5_hT3"  : 0.8100760357651001,
            "a^6_hT3"  : 0.3811852402868251,
            "a^7_hT3"  : 0.16419279860055377,
            "a^8_hT3"  : 0.06727158745187298,
            "a^9_hT3"  : 0.026738475289287966,
            "a^10_hT3" : 0.010426839830565283
        }
        self._params = [
            "a^0_hA1", "a^1_hA1", "a^2_hA1", "a^3_hA1", "a^4_hA1", "a^5_hA1", "a^6_hA1", "a^7_hA1", "a^8_hA1", "a^9_hA1", "a^10_hA1", 
            "a^0_hA2", "a^1_hA2", "a^2_hA2", "a^3_hA2", "a^4_hA2", "a^5_hA2", "a^6_hA2", "a^7_hA2", "a^8_hA2", "a^9_hA2", "a^10_hA2", 
            "a^0_hA3", "a^1_hA3", "a^2_hA3", "a^3_hA3", "a^4_hA3", "a^5_hA3", "a^6_hA3", "a^7_hA3", "a^8_hA3", "a^9_hA3", "a^10_hA3", 
            "a^0_hV", "a^1_hV", "a^2_hV", "a^3_hV", "a^4_hV", "a^5_hV", "a^6_hV", "a^7_hV", "a^8_hV", "a^9_hV", "a^10_hV", 
            "a^0_hT1", "a^1_hT1", "a^2_hT1", "a^3_hT1", "a^4_hT1", "a^5_hT1", "a^6_hT1", "a^7_hT1", "a^8_hT1", "a^9_hT1", "a^10_hT1", 
            "a^0_hT2", "a^1_hT2", "a^2_hT2", "a^3_hT2", "a^4_hT2", "a^5_hT2", "a^6_hT2", "a^7_hT2", "a^8_hT2", "a^9_hT2", "a^10_hT2", 
            "a^0_hT3", "a^1_hT3", "a^2_hT3", "a^3_hT3", "a^4_hT3", "a^5_hT3", "a^6_hT3", "a^7_hT3", "a^8_hT3", "a^9_hT3", "a^10_hT3"
        ]
        self._internalparams = {
            "Mb"  : self.par[f'm_{self.B}'],
            "Mc"  : self.par[f'm_{self.V}'],
            "lambdaqcdphys" : 0.5,
            "LambdaChi" : 1.0,
            "maxorder" : 10,
            "Mk" : 0.493677,
            "MPi" : 0.13957,
            "MEta" : 0.547862,
            "fpi" : 0.130
        }

    def _get_ff_h(self, q2: float, A: str) -> float:
        """Calculates h form factors according to https://arxiv.org/abs/2304.03137.
        Directly lifted from the ancillary files (LOAD_FIT.py)

        Parameters
        ----------
        q2 : float
        A : str
            FF to compute, hA1, hA2, hA3, hV, hT1, hT2, hT3

        Returns
        -------
        float
            FF computed at q2 value
        """
        value = 0
        w = self.w(q2)
        for order in range(int(self.internalparams["maxorder"]+1)):
            param_name = 'a^'+str(order)+'_'+A
            if param_name not in self.ffpar:
                print(f"{self.name} FF parameter {param_name} not found, not contributing")
                continue

            cumulator = self.ffpar[param_name]
            if order==0:
                value += cumulator
            else:
                value += cumulator*((w-1)**order)
        return value
    
    def get_ff(self, q2: float) -> dict[float]:
        """Implements FF caclculation according to https://arxiv.org/abs/2304.03137.
        Directly lifted from the ancillary files (LOAD_FIT.py)

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict[float]
            FF dictionary
        """
        mB = self.internalparams["Mb"]
        mV = self.internalparams["Mc"]
        r=mV/mB
        w = self.w(q2)
        hA1b=self._get_ff_h(q2, 'hA1')
        hA2b=self._get_ff_h(q2,'hA2')
        hA3b=self._get_ff_h(q2,'hA3')
        hVb =self._get_ff_h(q2,'hV')
        
        hT1b=self._get_ff_h(q2,'hT1')
        hT2b=self._get_ff_h(q2,'hT2')
        hT3b=self._get_ff_h(q2,'hT3')

        A0 = ((1.0+w)*hA1b+(r*w-1.0)*hA2b+(r-w)*hA3b)/(2.0*sqrt(r))
        A1 = hA1b*(1.0+w)*sqrt(r)/(1.0+r)
        A2 = (r*hA2b + hA3b)*(1.0+r)/(2.0*sqrt(r))
        V = hVb*(1.0+r)/(2.0*sqrt(r))

        T1b = -(1/(2*sqrt(r)))*(  (1-r)*hT2b - (1+r)*hT1b  )
        T2b = (1/(2*sqrt(r)))*(  2*r*(w+1)*hT1b/(1+r) - 2*r*(w-1)*hT2b/(1-r)  )
        T3b = (1/(2*sqrt(r)))*(  (1-r)*hT1b - (1+r)*hT2b + (1-r**2)*hT3b  )

        lambdabdstar = mB**4+mV**4+q2**2-2*(mB**2*mV**2+mB**2*q2+mV**2*q2)
        T23b = (1.0/(8*mB*mV**2))*((mB**2+3*mV**2-q2)*(mB+mV)*T2b - lambdabdstar*T3b/(mB-mV))

        ff = {
            "A0"  : A0,
            "A1"  : A1,
            "A2"  : A2,
            "V"   : V,
            "T1"  : T1b,
            "T2"  : T2b,
            "T3"  : T3b,
            "T23" : T23b,
        }
        ff['A12'] = ((ff['A1'] * (mB + mV)**2 * (mB**2 - mV**2 - q2)
                 - ff['A2'] * (mB**4 + (mV**2 - q2)**2
                 - 2 * mB**2 * (mV**2 + q2)))
                 / (16. * mB * mV**2 * (mB + mV)))

        return ff


