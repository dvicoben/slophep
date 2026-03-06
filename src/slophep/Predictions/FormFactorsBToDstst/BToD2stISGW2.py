import numpy as np
from slophep.Predictions.FormFactorsBToDstst import FormFactorBToD2st


class ISGW2_BToD2st(FormFactorBToD2st):
    def __init__(self, B: str, C: str, par: dict = None, scale: float = None, *ffargs) -> None:
        super().__init__(B, C, par, scale)
        self._name = "BToD2st_ISGW2"
        self._ffpar = {
        }
        self._params = []
        # from https://gitlab.com/mpapucci/Hammer/-/blob/v1.2.1/src/FormFactors/FFBtoD2starISGW2.cc
        internalparams = {
            "msb"     : 5.2                    ,
            "msd"     : 0.33                   ,
            "bb2"     : 0.431*0.431            ,
            "mbb"     : 5.31                   ,
            "msq"     : 1.82                   ,
            "bx2"     : 0.33*0.33              ,
            "mbx"     : (5.0*2.46+3.0*2.42)/8.0,
            "mqm"     : 0.1                    ,
            "nfp"     : 3.0                    ,
        }
        self._internalparams.update(internalparams)

        print(f"WARNING: {self.name} SM only, tensor FFs are zero.")
    
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

    def get_ff_mmeson(self, q2: float, mC: float, mB: float = None) -> dict:
        msb = self.internalparams["msb"]
        msd = self.internalparams["msd"]
        msq = self.internalparams["msq"]
        mbb = self.internalparams["mbb"]
        mbx = self.internalparams["mbx"]
        bb2 = self.internalparams["bb2"]
        bx2 = self.internalparams["bx2"]

        mb = self.internalparams["Mb"] if not mB else mB
        mx = mC
        mtb = msb + msd
        mtx = msq + msd
        mup = 1.0/(1.0/msq+1.0/msb)
        mum = 1.0/(1.0/msq-1.0/msb)
        bbx2 = 0.5*(bb2+bx2)
        tm = (mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        wt = 1.0+(tm-t)/(2.0*mbb*mbx)
        mqm = self.internalparams["mqm"]

        r2 = (
            3.0/(4.0*msb*msq)
            +3*msd*msd/(2*mbb*mbx*bbx2)
            +(16.0/(mbb*mbx*(33.0-2.0*self.internalparams["nfp"])))*np.log(self.get_as(mqm,mqm) / self.get_as(msq,msq))
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

    def get_ff(self, q2: float) -> dict:
        return self.get_ff_mmeson(q2, self.internalparams["Mc"], self.internalparams["Mb"])