#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 17:27:18 2019

@author: simon
"""

import tc_extract as tc
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
#sns.set_style('ticks')
params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','font.size':10,
          'xtick.labelsize':10, 'ytick.labelsize':10, 'axes.labelsize':10,
          'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True}
plt.rcParams.update(params)
from shapely.geometry import Point




###############################################################################
### This code will generate a pseudosection figure from a drawpd file
pseudo = tc.pseudo('dr-bg.txt')
f,a = pseudo.imitate_drawpd()


for linename in pseudo.tidyuni.keys():
    if pseudo.tidyuni[linename][1]==['liq']:
        a.plot(pseudo.tidyuni[linename][3],pseudo.tidyuni[linename][2],lw=2,c='k')
    if linename in pseudo.are[0][pseudo.tidyarea[-1]==['liq']]:
        a.plot(pseudo.tidyuni[linename][3],pseudo.tidyuni[linename][2],lw=2,c='k')

f.savefig('pseudo.pdf')

###############################################################################
### This code will run thermocalc over a grid

#gT, gP, phases = tc.makegrid(np.linspace(1100,2000,50),np.linspace(10.5,69.5,50))
#tcg = tcgridrun(gT,gP,phases)
#tcg.run()


###############################################################################
### This code will import the results from a grid run

# Toggle depending on whether its the first time importing data or not
#tcgr=tc.import_icfile('tc-kg1-ic.txt')
tcgr=tc.load_tcgrid()

#tcgr.results.to_csv('results.csv')
#results = pd.read_csv('results.csv')

if False:
    results = tcgr.results
    X = tcgr.gridresults()

    P = np.sort(results.pressure.unique())
    T = np.sort(results.temperature.unique())
    for x in range(np.shape(X)[0]):
        for y in range(np.shape(X)[1]):
            if pseudo.polygons[0][0].contains(Point(P[x],T[y])):
                X[x,y] = 1.0
            elif (X[x,y] > 0) == False:
                X[x,y] = 0.0

    f,a = plt.subplots()

    for linename in pseudo.tidyuni.keys():
        a.plot(pseudo.tidyuni[linename][3],pseudo.tidyuni[linename][2],lw=0.8,c='k',zorder=2)

    a.contourf(T,P,X,levels=[-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.999,1.0],cmap=plt.cm.Reds)

    plt.show()
