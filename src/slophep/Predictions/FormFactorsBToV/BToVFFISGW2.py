from math import sqrt
from slophep.Predictions.FormFactorsBToV import FormFactorBToV
import numpy as np


class ISGW2_BToV(FormFactorBToV):
    def __init__(self, B: str, V: str, par: dict = None, scale: float = None, *ffargs):
        super().__init__(B, V, par, scale)
        
        self._name = "BToV_ISGW2"
        self._ffpar = {
        }
        self._params = []

        internalparams = {
            "msb" : 5.2,
            "msd" : 0.33,
            "bb2" : 0.431*0.431,
            "mbb" : 5.31,
            "nf"  : 4.0,
            "cf"  : 0.989,
            "msq" : 1.82,
            "bx2" : 0.38*0.38,
            "mbx" : 0.75*2.01+0.25*1.87,
            "nfp" : 3.0,
            "mqm" : 0.1
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

    def get_ff(self, q2: float) -> dict:
        """ISGW2 FFs - in general not good

        Parameters
        ----------
        q2 : float
            q2 value to calculate FF at

        Returns
        -------
        dict
            FF dictionary
        """
        msb = self.internalparams["msb"]
        msd = self.internalparams["msd"]
        msq = self.internalparams["msq"]
        mbb = self.internalparams["mbb"]
        mbx = self.internalparams["mbx"]
        bb2 = self.internalparams["bb2"]
        bx2 = self.internalparams["bx2"]

        mtb=msb+msd
        mtx=msq+msd
        mup=1.0/(1.0/msq+1.0/msb)
        mum=1.0/(1.0/msq-1.0/msb)
        bbx2=0.5*(bb2+bx2)
        mb = self.internalparams["Mb"]
        mx = self.internalparams["Mc"]
        tm=(mb-mx)*(mb-mx)
        t = q2 if q2 < tm else 0.99*tm
        wt=1.0+(tm-t)/(2.0*mbb*mbx)
        mqm = self.internalparams["mqm"]

        r2=3.0/(4.0*msb*msq)+3*msd*msd/(2*mbb*mbx*bbx2) + (16.0/(mbb*mbx*(33.0-2.0*self.internalparams["nfp"])))*np.log(self.get_as(mqm,mqm) / self.get_as(msq,msq))
        ai = -1.0* ( 6.0/( 33.0 - 2.0*self.internalparams["nf"]))
        cji = np.power((self.get_as(msb,msb) / self.get_as(msq,msq)), ai)
        zji = msq / msb
        gammaji = self.get_gammaji(zji)
        chiji = -1.0 - ( gammaji / ( 1- zji ))
        
        betaji_g = (2.0/3.0)+gammaji
        betaji_f = (-2.0/3.0)+gammaji
        betaji_appam = -1.0-chiji+(4.0/(3.0*(1.0-zji)))+(2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))
        betaji_apmam = (1.0/3.0)-chiji-(4.0/(3.0*(1.0-zji)))-(2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))+gammaji

        r_g = cji*(1+(betaji_g*self.get_as( msq,sqrt(mb*msq) )/(np.pi)))
        r_f = cji*(1+(betaji_f*self.get_as( msq,sqrt(mb*msq) )/(np.pi)))
        r_apmam = cji*(1+(betaji_apmam*self.get_as( msq,sqrt(mb*msq) )/(np.pi)))
        
        f3=sqrt(mtx/mtb)*pow(sqrt(bx2*bb2)/bbx2,1.5)/((1.0+r2*(tm-t)/12.0)*(1.0+r2*(tm-t)/12.0))
        f3f=sqrt(mbx*mbb/(mtx*mtb))*f3
        f3g=sqrt(mtx*mtb/(mbx*mbb))*f3
        f3appam=sqrt(mtb*mtb*mtb*mbx/(mbb*mbb*mbb*mtx))*f3
        f3apmam=sqrt(mtx*mtb/(mbx*mbb))*f3

        cf = self.internalparams["cf"]
        f=cf*mtb*(1+wt+msd*(wt-1)/(2*mup))*f3f*r_f
        g=0.5*(1/msq-msd*bb2/(2*mum*mtx*bbx2))*f3g*r_g

        appam=cji*(msd*bx2*(1-msd*bx2/(2*mtb*bbx2))/
                    ((1+wt)*msq*msb*bbx2)-
                    betaji_appam*self.get_as( msq,sqrt(msq*mb) )/
                    (mtb*np.pi))*f3appam
        apmam=-1.0*(mtb/msb-msd*bx2/(2*mup*bbx2)+wt*msd*mtb*bx2*
                    (1-msd*bx2/(2*mtb*bbx2))/((wt+1)*msq*msb*bbx2))*f3apmam*r_apmam/mtx
        ap=0.5*(appam+apmam)
        am=0.5*(appam-apmam)
        # MbcSqq = pow(mb + mx, 2.) - t
        # fs = -2 * mx * f / MbcSqq
        
        # We need to fiddle around with some conversions to go from hammer to flavio basis
        r = mx/mb
        A1 = f/(mb + mx)
        pre = 1 / 2 / sqrt(mb * mx)
        A1_factor = pre * ((mb + mx)**2 - t) / (mb + mx)
        hA1 = A1/A1_factor
        hA2 = -(mb/sqrt(r))*(am + ap)
        hA3 = sqrt(r)*mb*(am - ap)
        A0 = pre * (((mb + mx)**2 - t) / (2 * mx) * hA1
                    - (mb**2 - mx**2 + t) / (2 * mb) * hA2
                    - (mb**2 - mx**2 - t) / (2 * mx) * hA3)
        A2 = pre * (mb + mx) * (hA3 + mx / mb * hA2)
        A12 = ((A1 * (mb + mx)**2 * (mb**2 - mx**2 - t)
                - A2 * (mb**4 + (mx**2 - t)**2
                - 2 * mb**2 * (mx**2 + t)))
                / (16. * mb * mx**2 * (mb + mx)))

        ff = {
            "V"   : (mb + mx) * g,
            "A0"  : A0,
            "A1"  : f/(mb + mx),
            "A2"  : A2,
            "A12" : A12,
            "T1"  : 0.0,
            "T2"  : 0.0,
            "T3"  : 0.0,
            "T23" : 0.0
        }
        return ff