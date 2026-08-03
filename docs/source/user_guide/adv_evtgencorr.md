# EvtGen Correction

As described in [arXiv:2602.18378](https://arxiv.org/abs/2602.18378), $B \to D^{**}\ell\nu$ modes generated in EvtGen suffer from artefacts that distort the phase-space distributions. In [arXiv:2602.18378](https://arxiv.org/abs/2602.18378), a corrective weighting is derived, which is a function of the $D^{**}$ mass ($m_R$) and has the form (shown for a 2-body $D^{**} \to h_1 h_2$, three-body modes have additional phase-space factors),

$$w(m_R) = m_R \frac{\Gamma_\mathrm{EvtGen}}{\Gamma} \int \mathrm{d}q^2 \frac{\mathrm{d}\Gamma}{\mathrm{d}q^2}\Bigg|_{m_R} $$

This is implemented in a small package by the authors of [arXiv:2602.18378](https://arxiv.org/abs/2602.18378), [mcambulance](https://github.com/Herren/MCAmbulance/tree/main) for a handful of decay modes and form-factor schemes. If the utilities in [mcambulance](https://github.com/Herren/MCAmbulance/tree/main) suit your purpose it is recommended you use that instead of SLOP's. Note that masses and widths in [mcambulance](https://github.com/Herren/MCAmbulance/tree/main) are taken from [`basf2`](https://github.com/belle2/basf2)'s [`evt.pdl`](https://github.com/belle2/basf2/blob/main/framework/particledb/data/evt.pdl) and this may not be the same in other EvtGen versions.

The implementation in SLOP produces an unnormalised correction, i.e. only the $m_R$ dependent portion, so by default SLOP's correction functions have the form (again, for a 2-body $D^{**} \to h_1 h_2$)

$$w_\mathrm{SLOP}(m_R) = m_R \int \mathrm{d}q^2 \frac{\mathrm{d}\Gamma}{\mathrm{d}q^2}\Bigg|_{m_R}$$

and it is up to the user to rescale the weights to the effective statistical size they desire. Note that $m_R$ is assumed to be in GeV in SLOP. An example of making this weighting function can be found [here](https://github.com/dvicoben/slophep/blob/master/python/evtgen_weighting.py).

For the phase-space factors in three-body modes, additional parameters are required (e.g. the width), and to avoid hard-coding them these have to be supplied by the user according to the contents of the `evt.pdl` in the version of EvtGen that was used to generate the sample. An example is also shown [here](https://github.com/dvicoben/slophep/blob/master/python/evtgen_weighting.py).