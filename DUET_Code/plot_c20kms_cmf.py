#!/usr/bin/env python3
import numpy as np
import seaborn as sb
from scipy.optimize import curve_fit
import plfit

import matplotlib.pyplot as plt

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetex=True)
#rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
rcParams.update({'font.size': 16})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True

masses = np.loadtxt('/Users/liuxihe/Desktop/我的/天协/科创数据/dendrogram_band6/sgrb1off_dendro_cores_band6.dat', usecols=(5), unpack=True)

fig, ax = plt.subplots(figsize=(8, 6))

logmasses = np.log10(masses)
print('min mass (Msun):',np.min(masses))
print('max mass (Msun):',np.max(masses))

sb.set_style('ticks')
cmfplot = sb.distplot(logmasses, hist=True, bins=20, rug=True, kde=False, color='k', hist_kws={'log':True}, rug_kws={'color':'r'}, norm_hist=False, ax=ax)
ax.set_xlim(np.log10([0.2,400]))
ax.set_ylim(0.9,80)
cmfplot.set_xticks(np.log10([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400]))
cmfplot.set_xticklabels(['','','','','','','','1','','','','','','','','','10','','','','','','','','','100','','',''])
cmfplot.set_yticks([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70])
cmfplot.set_yticklabels(['1','','','','','','','','','10','','','','','',''])
cmfplot.set_xlabel(r'Condensation Masses ($M_\odot$)')
cmfplot.set_ylabel(r'd$N$/d log $M$')
# Get the position for the label
x0, x1 = cmfplot.axes.get_xlim()
y0, y1 = cmfplot.axes.get_ylim()
labelx = x0 + (x1-x0)*0.58
labely = 10**(np.log10(y0) + np.log10(y1/y0)*0.85)
cmfplot.text(labelx,labely,r'CMZ Cloud ‘e’',size='large')
for vbar in np.log10([1,10,100]):
    cmfplot.axvline(vbar, color='gold', linewidth=1)

## Get data points and fit CMF
myfit = plfit.plfit(masses)

mlim = [myfit._xmin,np.max(masses)]
# The definition of alpha in the plfit package is different
index = myfit._alpha - 1.
error = myfit._alphaerr
#normal = 150.
n_core = myfit._nunique
normal = n_core / np.log10(mlim[1]-mlim[0])

def func_pl(x,a,b):
    return b*x**a

counts_fit = [func_pl(mlim[0],-1.*index,normal),func_pl(mlim[1],-1.*index,normal)]
#print(counts_fit)
tmp1 = str("%0.2f" % index)
tmp2 = str("%0.2f" % error)
plt.plot(np.log10(mlim),counts_fit,'m--',label=r'$\alpha=$ '+tmp1+r'$\pm$'+tmp2)
plt.legend(frameon=False, loc=(0.56,0.75))

# Mark the lower bound
x0, x1 = cmfplot.axes.get_xlim()
y0, y1 = cmfplot.axes.get_ylim()
labelx = x0 + (x1-x0)*0.626
labely = 10**(np.log10(y0) + np.log10(y1/y0)*0.70)
xmin = str('%0.2f' % myfit._xmin)
cmfplot.text(labelx,labely,r'$M_{min}=$ '+xmin+' $M_{\odot}$')

## Plot the mass sensitivity
cmfplot.axvline(np.log10(0.27), color='b', linestyle='--', linewidth=2)
cmfplot.axvline(np.log10(mlim[0]), color='magenta', linestyle='dotted', linewidth=2)

cmfplot.figure.savefig('core_mass.pdf',bbox_inches='tight')
