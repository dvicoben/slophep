import numpy as np
import flavio

# gammainfo = {
#     "D0*+" : 221.0e-3,
#     "D0*-" : 221.0e-3,
#     "D0*0" : 274.0e-3,
#     "D1*+" : 314.0e-3,
#     "D1*-" : 314.0e-3,
#     "D1*0" : 314.0e-3,
# }

def kallen_lambda(x: float, y: float, z: float) -> float:
    return x ** 2 + y ** 2 + z ** 2 - 2. * (x * y + y * z + z*x)


def BlattWeisskopf2(p: float, pp: float, R: float, l: float) -> float:
    if l == 0:
        return 1.
    elif l == 1:
        zR2 = R ** 2 * p ** 2
        zR2p = R ** 2 * pp ** 2
        return ( zR2p + 1. ) / ( zR2 + 1. )
    else:
        return 1.


def threeBodyPS(M2: float, M: float, m: float) -> float:
    # From MCAmbulance
    delta = (1. - (M + 2. * m) ** 2 / M2) / 8. / m / M
    pref = 32. * np.pi / np.sqrt(M / (M + 2. * m)) * (m * M) ** 3
    coefflist = np.array([0,0,1,m**2 + 8*m*M + M**2,(-5*m**4 + 100*m**3*M + 510*m**2*M**2 + 92*m*M**3 - M**4)/8., \
     (7*m**6 - 84*m**5*M + 987*m**4*M**2 + 4064*m**3*M**3 + 861*m**2*M**4 - 12*m*M**5 + M**6)/8., \
     (-105*m**8 + 1176*m**7*M - 7980*m**6*M**2 + 71624*m**5*M**3 + 259210*m**4*M**4 + 60200*m**3*M**5 - 876*m**2*M**6 + 120*m*M**7 - 9*M**8)/64., \
     (231*m**10 - 2640*m**9*M + 16335*m**8*M**2 - 82320*m**7*M**3 + 623550*m**6*M**4 + 2067696*m**5*M**5 + 509830*m**4*M**6 - 7280*m**3*M**7 + 1275*m**2*M**8 - \
      160*m*M**9 + 11*M**10)/64.,(-9009*m**12 + 108108*m**11*M - 671814*m**10*M**2 + 3031380*m**9*M**3 - 12615075*m**8*M**4 + 84849336*m**7*M**5 + \
      264017628*m**6*M**6 + 67908456*m**5*M**7 - 928935*m**4*M**8 + 189500*m**3*M**9 - 30774*m**2*M**10 + 3588*m*M**11 - 229*M**12)/1024., \
     (23595*m**14 - 300300*m**13*M + 1936935*m**12*M**2 - 8664656*m**11*M**3 + 31814783*m**10*M**4 - 115535420*m**9*M**5 + 711096155*m**8*M**6 + 2107726720*m**7*M**7 + \
      559484849*m**6*M**8 - 7258692*m**5*M**9 + 1646085*m**4*M**10 - 314160*m**3*M**11 + 47845*m**2*M**12 - 5236*m*M**13 + 313*M**14)/1024., \
     (-1042899*m**16 + 14119248*m**15*M - 95742504*m**14*M**2 + 440231792*m**13*M**3 - 1585315732*m**12*M**4 + 5018469456*m**11*M**5 - 16440191768*m**10*M**6 + \
      94415292400*m**9*M**7 + 269302755150*m**8*M**8 + 73253524080*m**7*M**9 - 897871128*m**6*M**10 + 220263120*m**5*M**11 - 47055764*m**4*M**12 + 8480752*m**3*M**13 - \
      1221672*m**2*M**14 + 126416*m*M**15 - 7123*M**16)/16384.,(3002285*m**18 - 43232904*m**17*M + 309697245*m**16*M**2 - 1485777696*m**15*M**3 + \
      5454036900*m**14*M**4 - 16778885760*m**13*M**5 + 47429549076*m**12*M**6 - 143253021600*m**11*M**7 + 778051467414*m**10*M**8 + 2151041754160*m**9*M**9 + \
      596651584230*m**8*M**10 - 6902314848*m**7*M**11 + 1799947380*m**6*M**12 - 418283712*m**5*M**13 + 84842820*m**4*M**14 - 14546400*m**3*M**15 + 1994661*m**2*M**16 - \
      196200*m*M**17 + 10469*M**18)/16384.])

    deltap = np.array([delta**k for k in range(len(coefflist))])
    res = np.dot(coefflist, deltap)
    res *= pref
    
    return res


def threeBody_LS2(M: float, 
                  M_nom: float, 
                  L: float, 
                  width: float,
                  daughters: list[str], 
                  par: dict[str, float]) -> tuple[float, float]:
    """
    Assumes 3 body is of form D** -> D pi pi
    """
    M2 = M*M
    gam = width
    masses = [par[f"m_{ip}"] for ip in daughters]
    m_D = np.max(masses)
    # Masses of non-D children
    m_other = [par[f"m_{ip}"] for ip in daughters if "D" not in ip]
    m = np.min(m_other)

    p  = threeBodyPS(M2, m_D, m)
    pp = 1.
    momratio = (p/pp) ** (2*L + 1)
    RBW = 3.

    ls2 = momratio / ((M2 - M_nom**2)**2 + M_nom**2 * gam**2) * BlattWeisskopf2(p, pp, RBW, L)
    ls2_nr = gam / 2. / np.pi / ( (np.sqrt(M2) - M_nom) ** 2 + gam ** 2 / 4. )
    return ls2, ls2_nr



def decrate_int(mres: float, obs) -> float:
    mB = obs.par['m_'+obs.B]
    q2min = obs.q2min
    q2max = (mB - mres)**2
    def integrand(qsq):
        rate = obs.dGdq2M(qsq, mres)
        return rate
    return flavio.math.integrate.nintegrate(integrand, q2min, q2max)

