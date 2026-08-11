# Copyright (C) 2026  David Vico Benet

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# SLOP or SLOPHEP employs, translates and/or reimplements utilities from:
# - flavio (https://flav-io.github.io/), which is distributed under the MIT License, 
# and without any warranty, see <https://mit-license.org/>
# - Hammer (https://hammer.physics.lbl.gov/), which is distributed under version 3 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>
# - EOS (https://eoshep.org/), which is distributed under version 2 of the GPL, 
# and without any warranty, see <https://www.gnu.org/licenses/>

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