"""
Bs->Ds* Form-factors
"""
from typing import Any
import numpy as np

import slophep.FormFactors.FormFactorsBToV as FFBToV
from slophep.Tools.errfluct_tools import fluctsettings, FluctType
from slophep.Core.user_registry import FFregistry

import logging
logger = logging.getLogger(__name__)

@FFregistry.register
class ISGW2(FFBToV.FFBToV_ISGW2):
    _name = "FFBsToDsst@ISGW2"
    def __init__(self):
        super().__init__("Bs", "Ds*")

    def define_userparams(self) -> dict[str, Any]:
        ffpar = { 
            "msb" : 5.2,
            "msd" : 0.55,
            "bb2" : 0.54*0.54,
            "mbb" : 5.38,
            "nf"  : 4.0,
            "cf"  : 0.984,
            "msq" : 1.82,
            "bx2" : 0.49*0.49,
            "mbx" : 0.75*2.11+0.25*1.97,
            "nfp" : 3.0,
            "mqm" : 0.1
        }
        return ffpar
        

@FFregistry.register
class BLPR(FFBToV.FFBToV_BLPR):
    _name = "FFBsToDsst@BLPR"
    def __init__(self):
        super().__init__("Bs", "Ds*")


@FFregistry.register
class BLPRXP(FFBToV.FFBToV_BLPRXP):
    _name = "FFBsToDsst@BLPRXP"
    def __init__(self):
        super().__init__("Bs", "Ds*")
    

@FFregistry.register
class CLN(FFBToV.FFBToV_CLN):
    _name = "FFBsToDsst@CLN"
    def __init__(self):
        super().__init__("Bs", "Ds*")


@FFregistry.register
class CLN2(FFBToV.FFBToV_CLN2):
    _name = "FFBsToDsst@CLN2"
    def __init__(self):
        super().__init__("Bs", "Ds*")


@FFregistry.register
class BGL(FFBToV.FFBToV_BGL):
    _name = "FFBsToDsst@BGL"
    def __init__(self):
        logger.info("This is the same as BdToDst but with different meson masses")
        super().__init__("Bs", "Ds*")


@FFregistry.register
class HPQCD2023(FFBToV.FFBToV_HPQCD):
    _name = "FFBsToDsst@HPQCD2023"
    def __init__(self):
        super().__init__("Bs", "Ds*")

    def define_userparams(self) -> dict[str, Any]:
        pars = super().define_userparams()
        pars.update({
            "s_a^0_hA1_Mpi^2/Lambda" : -0.03126051843850383,
            "s_a^1_hA1_Mpi^2/Lambda" : 0.3395753129452487,
            "s_a^2_hA1_Mpi^2/Lambda" : -0.0378096661736465,
            "s_a^3_hA1_Mpi^2/Lambda" : -0.00022588027772022152,
            "s_a^4_hA1_Mpi^2/Lambda" : -0.0009771511385570624,
            "s_a^5_hA1_Mpi^2/Lambda" : 0.0009438671414110913,
            "s_a^6_hA1_Mpi^2/Lambda" : -0.0004853084408438801,
            "s_a^7_hA1_Mpi^2/Lambda" : 0.000195564469102004,
            "s_a^8_hA1_Mpi^2/Lambda" : -6.860628675070964e-05,
            "s_a^9_hA1_Mpi^2/Lambda" : 2.1958839264703956e-05,
            "s_a^10_hA1_Mpi^2/Lambda": -6.585624503198971e-06,
            "s_a^0_hA2_Mpi^2/Lambda" : 0.3182874441028258,
            "s_a^1_hA2_Mpi^2/Lambda" : -0.002358007805452479,
            "s_a^2_hA2_Mpi^2/Lambda" : -0.036131937656837684,
            "s_a^3_hA2_Mpi^2/Lambda" : 0.0200461048107326,
            "s_a^4_hA2_Mpi^2/Lambda" : -0.007319302649209761,
            "s_a^5_hA2_Mpi^2/Lambda" : 0.0022640881937483863,
            "s_a^6_hA2_Mpi^2/Lambda" : -0.0006462286490616317,
            "s_a^7_hA2_Mpi^2/Lambda" : 0.00017688858062101217,
            "s_a^8_hA2_Mpi^2/Lambda" : -4.734916184639905e-05,
            "s_a^9_hA2_Mpi^2/Lambda" : 1.2517774155015858e-05,
            "s_a^10_hA2_Mpi^2/Lambda": -3.283828236776756e-06,
            "s_a^0_hA3_Mpi^2/Lambda" : -0.4669421143513379,
            "s_a^1_hA3_Mpi^2/Lambda" : -0.3007476956387983,
            "s_a^2_hA3_Mpi^2/Lambda" : 0.04611470661424419,
            "s_a^3_hA3_Mpi^2/Lambda" : 0.011521010278445149,
            "s_a^4_hA3_Mpi^2/Lambda" : -0.008583425778290835,
            "s_a^5_hA3_Mpi^2/Lambda" : 0.0033026661287221745,
            "s_a^6_hA3_Mpi^2/Lambda" : -0.0010446249983807135,
            "s_a^7_hA3_Mpi^2/Lambda" : 0.00030191005903831754,
            "s_a^8_hA3_Mpi^2/Lambda" : -8.314760253633629e-05,
            "s_a^9_hA3_Mpi^2/Lambda" : 2.225556459295046e-05,
            "s_a^10_hA3_Mpi^2/Lambda": -5.845426029712506e-06,
            "s_a^0_hV_Mpi^2/Lambda"  : 0.12471461687559403,
            "s_a^1_hV_Mpi^2/Lambda"  : -0.05004397957696898,
            "s_a^2_hV_Mpi^2/Lambda"  : 0.1289606938018353,
            "s_a^3_hV_Mpi^2/Lambda"  : -0.079086309038104,
            "s_a^4_hV_Mpi^2/Lambda"  : 0.032844043697673075,
            "s_a^5_hV_Mpi^2/Lambda"  : -0.011562669642619974,
            "s_a^6_hV_Mpi^2/Lambda"  : 0.0037251139092783483,
            "s_a^7_hV_Mpi^2/Lambda"  : -0.0011362451270477342,
            "s_a^8_hV_Mpi^2/Lambda"  : 0.000334071638565463,
            "s_a^9_hV_Mpi^2/Lambda"  : -9.568986504251624e-05,
            "s_a^10_hV_Mpi^2/Lambda" : 2.6887324046381517e-05,
            "s_a^0_hT1_Mpi^2/Lambda" : 0.03664848541553486,
            "s_a^1_hT1_Mpi^2/Lambda" : 0.18657682691589345,
            "s_a^2_hT1_Mpi^2/Lambda" : -0.08093365619055709,
            "s_a^3_hT1_Mpi^2/Lambda" : -0.030112048327859926,
            "s_a^4_hT1_Mpi^2/Lambda" : 0.008278305294355453,
            "s_a^5_hT1_Mpi^2/Lambda" : 0.002340592610290264,
            "s_a^6_hT1_Mpi^2/Lambda" : 0.0003882344314551288,
            "s_a^7_hT1_Mpi^2/Lambda" : 4.7734956111884964e-05,
            "s_a^8_hT1_Mpi^2/Lambda" : 4.2104640226983955e-06,
            "s_a^9_hT1_Mpi^2/Lambda" : 1.0420513102939777e-07,
            "s_a^10_hT1_Mpi^2/Lambda": -6.49929920573398e-08,
            "s_a^0_hT2_Mpi^2/Lambda" : 0.0952209853630684,
            "s_a^1_hT2_Mpi^2/Lambda" : 0.001613489830598449,
            "s_a^2_hT2_Mpi^2/Lambda" : -0.01463914475311275,
            "s_a^3_hT2_Mpi^2/Lambda" : 0.11116434089187299,
            "s_a^4_hT2_Mpi^2/Lambda" : 0.022324639963270988,
            "s_a^5_hT2_Mpi^2/Lambda" : 0.0034425833737695773,
            "s_a^6_hT2_Mpi^2/Lambda" : 0.0004757339143352295,
            "s_a^7_hT2_Mpi^2/Lambda" : 6.26028538461786e-05,
            "s_a^8_hT2_Mpi^2/Lambda" : 8.070907701985311e-06,
            "s_a^9_hT2_Mpi^2/Lambda" : 1.0334726477604513e-06,
            "s_a^10_hT2_Mpi^2/Lambda": 1.322628801385152e-07,
            "s_a^0_hT3_Mpi^2/Lambda" : 0.006087770769255353,
            "s_a^1_hT3_Mpi^2/Lambda" : -0.003161752641004126,
            "s_a^2_hT3_Mpi^2/Lambda" : -0.0013687749514931537,
            "s_a^3_hT3_Mpi^2/Lambda" : 0.0026863645489903377,
            "s_a^4_hT3_Mpi^2/Lambda" : 0.0013330243250929425,
            "s_a^5_hT3_Mpi^2/Lambda" : 0.00033830691934960124,
            "s_a^6_hT3_Mpi^2/Lambda" : 6.560748720422526e-05,
            "s_a^7_hT3_Mpi^2/Lambda" : 1.1088161729243377e-05,
            "s_a^8_hT3_Mpi^2/Lambda" : 1.7326179440838514e-06,
            "s_a^9_hT3_Mpi^2/Lambda" : 2.581703455837798e-07,
            "s_a^10_hT3_Mpi^2/Lambda": 3.733799623256207e-08,
            #internalparams
            "LambdaChi" : 1.0,
            "Mk"        : 0.493677,
            "MPi"       : 0.13957,
            "MEta"      : 0.547862,
            "fpi"       : 0.130
        })
        return pars

    def _get_ff_h(self, q2: float, A: str) -> float:
            value = 0
            w = self.w(q2)
            mK = self.get_userparam("Mk")
            mPI = self.get_userparam("MPi")
            LambdaChi = self.get_userparam("LambdaChi")
            chi_mult = ((mK/LambdaChi)**2-(mPI/LambdaChi)**2)
            for order in range(int(self.get_userparam("maxorder")+1)):
                param_name = 'a^'+str(order)+'_'+A
                if f"{self.name}:{param_name}" not in self.pm.params:
                    logger.info(f"{self.name} FF parameter {param_name} not found, not contributing")
                    continue
    
                cumulator = self.get_userparam(param_name)
                chicumulator = self.get_userparam(f"s_{param_name}_Mpi^2/Lambda")*chi_mult
                if order==0:
                    value += cumulator*(1+chicumulator)
                else:
                    value += cumulator*((w-1)**order)*(1+chicumulator)
            return value



@FFregistry.register
class BGL_Hammer(FFBToV.FFBToV_BGLGeneric):
    _name = "FFBsToDsst@BGL_Hammer"
    def __init__(self):
        super().__init__("Bs", "Ds*", 3, 3, 3, 2, 0, 0, 0)

    def define_userparams(self) -> dict[str, Any]:
        # Meant to reproduce https://arxiv.org/pdf/2105.14019
        # NOTE: c0 is fixed by kinematical constraint in Eq. (72)
        pars = super().define_userparams()
        ffpar = {
            "a_g_0"  : 0.00038   ,
            "a_g_1"  : 0.026905  ,
            "a_g_2"  : 0.        ,
            "a_f_0"  : 0.00055   ,
            "a_f_1"  : -0.0020370,
            "a_f_2"  : 0.        ,
            "a_F1_1" : -0.000433 ,
            "a_F1_2" : 0.005353  ,
            "a_F2_0" : 0.007     ,
            "a_F2_1" : -0.036    ,
            # internalparams
            "Vcb"        : 41.5e-3 ,                       
            "chiF1"      : 3.068e-4, # GeV^-2
            "chif"       : 3.068e-4, # GeV^-2
            "chig"       : 5.280e-4, # GeV^-2
            "chiF2"      : 2.466e-3,
            "nc"         : 2.6     ,
            "etaEW"      : 1.0066  ,
            "BcStates1p" : np.array([6.730, 6.736, 7.135, 7.142]),  # GeV
            "BcStates1m" : np.array([6.337, 6.899, 7.012, 7.280]),  # GeV
            "BcStates0m" : np.array([6.275, 6.842, 7.250])          # GeV
        }
        pars.update(ffpar)
        return pars

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff(self, q2: float) -> dict:
        """Calculates BGL form factors (SM only) as in hammer https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoDstarBGL.cc?ref_type=tags
        
        Note that there is an additional 1./(etaEW*Vcb) factor applied to FFs as in Hammer.

        Parameters
        ----------
        q2 : float

        Returns
        -------
        dict
            FF dictionary
        """
        etaEWVcb = self.get_userparam("etaEW")*self.get_userparam("Vcb")
        ff = super().get_ff(q2)
        return {k : ff[k]/etaEWVcb for k in ff}


@FFregistry.register
class BLPR_Hammer(FFBToV.FFBToV_BLPR):
    _name = "FFBsToDsst@BLPR_Hammer"
    def __init__(self):
        super().__init__("Bs", "Ds*")

    def define_userparams(self):
        ffpar = {
            "RhoSq"     : 1.24      ,
            "Chi21"     : -0.06     ,
            "Chi2p"     : 0.0       ,
            "Chi3p"     : 0.05      ,
            "Eta1"      : 0.30      ,
            "Etap"      : -0.05     ,
            "dV20"      : 0.0       ,
            # internalparams:
            "mb"        : 4.710     ,
            "delta_mbc" : 3.4       ,
            "normscale" : 1.0       ,
            "ash"       : 0.26/np.pi,
            "ebReb"     : 0.861     ,
            "ecRec"     : 0.822     ,
            "mbarB"     : 5.313     ,
            "lam1"      : -0.3      , # GeV^2
        }
        return ffpar


@FFregistry.register
class CLN_Hammer(FFBToV.FFBToV_CLN):
    _name = "FFBsToDsst@CLN_Hammer"
    def __init__(self):
        super().__init__("Bs", "Ds*")

    def define_userparams(self):
        ffpar = {
            "RhoSq" : 1.207,
            "h_A1"  : 0.908,
            "R1"    : 1.401,
            "R2"    : 0.854,
            "R0"    : 1.15
        }
        return ffpar