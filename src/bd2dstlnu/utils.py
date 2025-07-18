import matplotlib.pyplot as plt

def setPlotParams(params = {}):
    """Set of parameters for plots in matplotlib
    """
    if len(params) > 0:
        plt.rcParams.update(params)
    else:
        font = {'family' : 'cmr10', #Or Times New Roman
            'size'   : 14}
        plt.rc('font', **font)
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams["axes.formatter.use_mathtext"] = True
        plt.rcParams.update({'axes.labelsize': 14,
                            'legend.fontsize': 10,
                            'xtick.labelsize': 12,
                            'ytick.labelsize': 12,
                            'figure.figsize': [8, 8/1.618]})