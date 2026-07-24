import numpy as np
from slophep.FormFactors.FormFactorsBToDstst.FFBToDststBase import (
    FormFactorBToD0st, FormFactorBToD1, FormFactorBToD1st, FormFactorBToD2st
)
from slophep.Tools.SamplingTools import fluctsettings, FluctType

import logging
logger = logging.getLogger(__name__)


class FFBToD1_ISGW2(FormFactorBToD1):
    _name = "FFBToD1@ISGW2"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "msb"     : 5.2                    ,
            "msd"     : 0.33                   ,
            "bb2"     : 0.431*0.431            ,
            "mbb"     : 5.31                   ,
            "msq"     : 1.82                   ,
            "bx2"     : 0.33*0.33              ,
            "mbx"     : (5.0*2.46+3.0*2.42)/8.0,
            "mqm"     : 0.1                    ,
            "nfp"     : 3.0                    ,
            "SmearQ2" : True
        }
        return ffpar

    def get_gammaji(self, z: float) -> float:
            value = 2+((2.0*z)/(1-z))*np.log(z)
            return -1.0*value
    
    def get_as(self, mq1: float, mq2: float) -> float:
        if mq2 <= 0.6:
            return 0.6
        lambdaSq = 0.04
        Nf = 4 if mq1 >= 1.85 else 3.0
        value = 12.0*np.pi/(33.0-2.0*Nf)
        value /= np.log(mq2*mq2/lambdaSq)
        return value

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        msb = self.get_userparam("msb")
        msd = self.get_userparam("msd")
        msq = self.get_userparam("msq")
        mbb = self.get_userparam("mbb")
        mbx = self.get_userparam("mbx")
        bb2 = self.get_userparam("bb2")
        bx2 = self.get_userparam("bx2")

        mtb = msb + msd
        mtx = msq + msd
        mup = 1.0/(1.0/msq+1.0/msb)
        mum = 1.0/(1.0/msq-1.0/msb)
        bbx2 = 0.5*(bb2+bx2)
        mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        mx = mC
        tm = (mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        wt = 1.0+(tm-t)/(2.0*mbb*mbx)
        mqm = self.get_userparam("mqm")
        nfp = self.get_userparam("nfp")

        r2 = (
            3.0/(4.0*msb*msq)
            +3*msd*msd/(2*mbb*mbx*bbx2)
            +(16.0/(mbb*mbx*(33.0-2.0*nfp)))*np.log(self.get_as(mqm,mqm) / self.get_as(msq,msq))
        )
        f5 = np.sqrt(mtx/mtb)*np.power(np.sqrt(bx2*bb2)/bbx2,5.0/2.0)/(np.power((1.0+r2*(tm-t)/18.0), 3.0))
        f5v = f5*np.power(( mbb / mtb ),-0.5)*np.power((mbx/mtx),-0.5)
        f5r = f5*np.power(( mbb / mtb ),0.5)*np.power((mbx/mtx),0.5)
        f5sppsm = f5*np.power(( mbb / mtb ),-1.5)*np.power((mbx/mtx),0.5)
        f5spmsm = f5*np.power(( mbb / mtb ),-0.5)*np.power((mbx/mtx),-0.5)

        vv = -msd*f5v/(2.0*np.sqrt(3.0*bb2)*mtx)*((wt+1)/2.0+bb2*mtb/(2.0*msd*msq*msb))
        rr = -2.*mtb*np.sqrt(bb2/3.0)*f5r*(1.0/msq + mtx*msd*(wt-1)/(2.0*bb2)*((wt+1)/(2.0*msq)-msd*bb2/(2.0*mum*mtx*bbx2)))
        sppsm = -np.sqrt(3.0)*msd*f5sppsm/(2.0*np.sqrt(bb2)*mtb)*(1 - msd/(3.0*msq) - msd*bb2/(3.0*bbx2)*(1.0/(2.0*mum)-1.0/mup))
        spmsm = -msd*f5spmsm/(2.0*np.sqrt(3.0*bb2)*mtx)*((2-wt)*mtx/msq + msd*bb2/bbx2*(1.0/(2.0*mum)-1.0/mup))

        sqmbmx = np.sqrt(mb*mx)
        Fv1 = rr/sqmbmx
        Fv2 = mb*mb*sppsm/sqmbmx
        Fv3 = spmsm*sqmbmx
        Fa = 2*vv*sqmbmx

        # Ad hoc smearing, following https://repo.hepforge.org/source/evtgen/browse/master/src/EvtGenModels/EvtISGW2FF.cpp
        smearQ2 = 1.0
        if self.get_userparam("SmearQ2"):
            q2max = (mb-mx)**2
            q2maxmean = (self.get_param(f"m_{self.B}")-self.get_userparam(f"m_{self.C}"))**2
            smearQ2 = np.min((np.sqrt(q2maxmean/q2max), 1000.))

        ffs = {
            "fS"  : 0.0,
            "fV1" : smearQ2*Fv1,
            "fV2" : smearQ2*Fv2,
            "fV3" : smearQ2*Fv3,
            "fA"  : smearQ2*Fa,
            "fT1" : 0.0,
            "fT2" : 0.0,
            "fT3" : 0.0
        }
        
        return ffs



class FFBToD2st_ISGW2(FormFactorBToD2st):
    _name = "FFBToD2st@ISGW2"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "msb" : 5.2                    ,
            "msd" : 0.33                   ,
            "bb2" : 0.431*0.431            ,
            "mbb" : 5.31                   ,
            "msq" : 1.82                   ,
            "bx2" : 0.33*0.33              ,
            "mbx" : (5.0*2.46+3.0*2.42)/8.0,
            "mqm" : 0.1                    ,
            "nfp" : 3.0                    ,
        }
        return ffpar

    def get_gammaji(self, z: float) -> float:
            value = 2+((2.0*z)/(1-z))*np.log(z)
            return -1.0*value
    
    def get_as(self, mq1: float, mq2: float) -> float:
        if mq2 <= 0.6:
            return 0.6
        lambdaSq = 0.04
        Nf = 4 if mq1 >= 1.85 else 3.0
        value = 12.0*np.pi/(33.0-2.0*Nf)
        value /= np.log(mq2*mq2/lambdaSq)
        return value

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        msb = self.get_userparam("msb")
        msd = self.get_userparam("msd")
        msq = self.get_userparam("msq")
        mbb = self.get_userparam("mbb")
        mbx = self.get_userparam("mbx")
        bb2 = self.get_userparam("bb2")
        bx2 = self.get_userparam("bx2")

        mtb = msb + msd
        mtx = msq + msd
        mup = 1.0/(1.0/msq+1.0/msb)
        mum = 1.0/(1.0/msq-1.0/msb)
        bbx2 = 0.5*(bb2+bx2)
        mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        mx = mC
        tm = (mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        wt = 1.0+(tm-t)/(2.0*mbb*mbx)
        mqm = self.get_userparam("mqm")
        nfp = self.get_userparam("nfp")

        r2 = (
            3.0/(4.0*msb*msq)
            +3*msd*msd/(2*mbb*mbx*bbx2)
            +(16.0/(mbb*mbx*(33.0-2.0*nfp)))*np.log(self.get_as(mqm,mqm) / self.get_as(msq,msq))
        )
        
        f5 = np.sqrt(mtx/mtb)*np.power(np.sqrt(bx2*bb2)/bbx2,5.0/2.0)/(np.power((1.0+r2*(tm-t)/18.0),3.0))

        f5h = f5*np.power(( mbb / mtb ),-1.5)*np.power((mbx/mtx),-0.5)
        f5k = f5*np.power(( mbb / mtb ),-0.5)*np.power((mbx/mtx),0.5)
        f5bppbm = f5*np.power(( mbb / mtb ),-2.5)*np.power((mbx/mtx),0.5)
        f5bpmbm = f5*np.power(( mbb / mtb ),-1.5)*np.power((mbx/mtx),-0.5)

        hf = f5h*(msd/(np.sqrt(8.0*bb2)*mtb))*((1.0/msq)-(msd*bb2/(2.0*mum*mtx*bbx2)))
        kf = f5k*(msd/(np.sqrt(2.0*bb2)))*(1.0+wt)
        bppbm = ((msd*msd*f5bppbm*bx2)/(np.sqrt(32.0*bb2)*msq*msb*mtb*bbx2))*(1.0-(msd*bx2/(2.0*mtb*bbx2)))
        bpmbm = -1.0*(msd*f5bpmbm/(np.sqrt(2.0*bb2)*msb*mtx))*(
            1.0 - ((msd*msb*bx2)/(2.0*mup*mtb*bbx2))
            + (
                (msd*bx2*(1.0-((msd*bx2)/(2.0*mtb*bbx2))))
                / (4.0*msq*bbx2)
            )
        )

        sqmbmx = np.sqrt(mb*mx)
        Ka1 = kf*mb/sqmbmx
        Ka2 = mb*mb*mb*bppbm/sqmbmx
        Ka3 = bpmbm*mb*sqmbmx
        Kv = 2*hf*mb*sqmbmx

        ffs = {
            "kP"  : 0.0,
            "kA1" : Ka1,
            "kA2" : Ka2,
            "kA3" : Ka3,
            "kV"  : Kv,
            "kT1" : 0.0,
            "kT2" : 0.0,
            "kT3" : 0.0,
        }
        return ffs



class FFBToD0st_ISGW2(FormFactorBToD0st):
    _name = "FFBToD0st@ISGW2"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "msb" : 5.2                    ,
            "msd" : 0.33                   ,
            "bb2" : 0.431*0.431            ,
            "mbb" : 5.31                   ,
            "msq" : 1.82                   ,
            "bx2" : 0.33*0.33              ,
            "mbx" : (3.0*2.49+2.40)/4.0    ,
            "mqm" : 0.1                    ,
            "nfp" : 3.0
        }
        return ffpar

    def get_gammaji(self, z: float) -> float:
            value = 2+((2.0*z)/(1-z))*np.log(z)
            return -1.0*value
    
    def get_as(self, mq1: float, mq2: float) -> float:
        if mq2 <= 0.6:
            return 0.6
        lambdaSq = 0.04
        Nf = 4 if mq1 >= 1.85 else 3.0
        value = 12.0*np.pi/(33.0-2.0*Nf)
        value /= np.log(mq2*mq2/lambdaSq)
        return value

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        msb = self.get_userparam("msb")
        msd = self.get_userparam("msd")
        msq = self.get_userparam("msq")
        mbb = self.get_userparam("mbb")
        mbx = self.get_userparam("mbx")
        bb2 = self.get_userparam("bb2")
        bx2 = self.get_userparam("bx2")

        mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        mx = mC
        mtb = msb+msd
        mtx = msq+msd
        bbx2 = 0.5*(bb2+bx2)
        
        tm = (mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        mqm = self.get_userparam("mqm")
        nfp = self.get_userparam("nfp")

        r2 = (
            3.0/(4.0*msb*msq)
            +3*msd*msd/(2*mbb*mbx*bbx2) 
            + (16.0/(mbb*mbx*(33.0-2.0*nfp)))*np.log(self.get_as(mqm, mqm)/self.get_as(msq, msq))
        )
        f5 = np.sqrt(mtx/mtb)*np.power(np.sqrt(bx2*bb2)/bbx2,5.0/2.0)/(np.power((1.0+r2*(tm-t)/18.0),3.0))
        
        f5uppum = f5*np.power((mbb/mtb), -0.5)*np.power((mbx/mtx), 0.5)
        f5upmum = f5*np.power((mbb/mtb), 0.5)*np.power((mbx/mtx), -0.5)

        uppum = -1.0*f5uppum*np.sqrt(2.0/(3.0*bb2))*msd
        upmum = 1.0*f5upmum*np.sqrt(2.0/(3.0*bb2))*msd*mtb/mtx

        gppgm = mb/np.sqrt(mb*mx)*uppum
        gpmgm = mx/np.sqrt(mb*mx)*upmum
        gp = (gppgm + gpmgm)/2.0
        gm = (gppgm - gpmgm)/2.0


        ffs = {
            "gP"  : 0.0,
            "g+"  : gp ,
            "g-"  : gm ,
            "gT"  : 0.0
        }
        return ffs



class FFBToD1st_ISGW2(FormFactorBToD1st):
    _name = "FFBToD1st@ISGW2"
    def __init__(self, B: str, C: str):
        logger.info(f"{self.name} tensor FFs are zero.")
        super().__init__(B, C)

    def define_userparams(self):
        ffpar = {
            "msb"     : 5.2                    ,
            "msd"     : 0.33                   ,
            "bb2"     : 0.431*0.431            ,
            "mbb"     : 5.31                   ,
            "msq"     : 1.82                   ,
            "bx2"     : 0.33*0.33              ,
            "mbx"     : (3.0*2.49+2.40)/4.0    ,
            "mqm"     : 0.1                    ,
            "nfp"     : 3.0                    ,
            "SmearQ2" : True
        }
        return ffpar

    def get_gammaji(self, z: float) -> float:
            value = 2+((2.0*z)/(1-z))*np.log(z)
            return -1.0*value
    
    def get_as(self, mq1: float, mq2: float) -> float:
        if mq2 <= 0.6:
            return 0.6
        lambdaSq = 0.04
        Nf = 4 if mq1 >= 1.85 else 3.0
        value = 12.0*np.pi/(33.0-2.0*Nf)
        value /= np.log(mq2*mq2/lambdaSq)
        return value

    @fluctsettings(FluctType.DICTNUMERIC)
    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        msb = self.get_userparam("msb")
        msd = self.get_userparam("msd")
        msq = self.get_userparam("msq")
        mbb = self.get_userparam("mbb")
        mbx = self.get_userparam("mbx")
        bb2 = self.get_userparam("bb2")
        bx2 = self.get_userparam("bx2")

        mb = mB if mB is not None else self.get_param(f"m_{self.B}")
        mx = mC
        mtb = msb+msd
        mtx = msq+msd
        mum = 1.0/(1.0/msq-1.0/msb)
        bbx2 = 0.5*(bb2+bx2)
        
        tm = (mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        wt = 1.0+(tm-t)/(2.0*mbb*mbx)
        mqm = self.get_userparam("mqm")
        nfp = self.get_userparam("nfp")

        r2 = (
            3.0/(4.0*msb*msq)
            +3*msd*msd/(2*mbb*mbx*bbx2) 
            + (16.0/(mbb*mbx*(33.0-2.0*self.internalparams["nfp"])))*np.log(self.get_as(mqm, mqm)/self.get_as(msq, msq))
        )
        f5 = np.sqrt(mtx/mtb)*np.power(np.sqrt(bx2*bb2)/bbx2,5.0/2.0)/(np.power((1.0+r2*(tm-t)/18.0),3.0))
        f5q = f5*np.power(( mbb / mtb ),-0.5)*np.power((mbx/mtx),-0.5)
        f5l = f5*np.power(( mbb / mtb ),0.5)*np.power((mbx/mtx),0.5)
        f5cppcm = f5*np.power(( mbb / mtb ),-1.5)*np.power((mbx/mtx),0.5)
        f5cpmcm = f5*np.power(( mbb / mtb ),-0.5)*np.power((mbx/mtx),-0.5)

        ql = f5q*np.sqrt(1.0/6.0)*msd/(np.sqrt(bb2)*mtx)*(1.0-bb2*mtb/(4.0*msd*msq*msb))
        ll = f5l*np.sqrt(2.0/3.0)*mtb*np.sqrt(bb2)*(1.0/(2.0*msq) - 3.0/(2.0*msb) + msd*mtx*(wt-1)/bb2*(1.0/msq-msd*bb2/(2.0*mum*mtx*bbx2)))
        cppcm = msd*msd*bx2*f5cppcm/(np.sqrt(6.0)*mtb*msq*np.sqrt(bb2)*bbx2)
        cpmcm = -np.sqrt(2.0/3.0)*msd*f5cpmcm/(np.sqrt(bb2)*mtx)*(1+msd*bx2/(2.0*msq*bbx2))

        sqmbmx = np.sqrt(mb*mx)
        Gv1 = ll/sqmbmx
        Gv2 = mb*mb*cppcm/sqmbmx
        Gv3 = cpmcm*sqmbmx
        Ga = 2*ql*sqmbmx

        # Ad hoc smearing, following https://repo.hepforge.org/source/evtgen/browse/master/src/EvtGenModels/EvtISGW2FF.cpp
        smearQ2 = 1.0
        if self.get_userparam("SmearQ2"):
            q2max = (mb-mx)**2
            q2maxmean = (self.get_param(f"m_{self.B}")-self.get_userparam(f"m_{self.C}"))**2
            smearQ2 = np.min((np.sqrt(q2maxmean/q2max), 1000.))

        ffs = {
            "gS"  : 0.0,
            "gV1" : smearQ2*Gv1,
            "gV2" : smearQ2*Gv2,
            "gV3" : smearQ2*Gv3,
            "gA"  : smearQ2*Ga,
            "gT1" : 0.0,
            "gT2" : 0.0,
            "gT3" : 0.0
        }
        return ffs



