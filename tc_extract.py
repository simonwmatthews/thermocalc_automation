#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 23:00:53 2019

@author: simon
"""

from scipy.spatial import Delaunay
from scipy.interpolate import interp1d
from scipy.optimize import newton
from scipy.spatial import ConvexHull
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry import Point
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
#sns.set_style('ticks')
params = {'text.usetex': False, 'mathtext.fontset': 'stixsans','font.size':10,
          'xtick.labelsize':10, 'ytick.labelsize':10, 'axes.labelsize':10,
          'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True}
plt.rcParams.update(params)
import subprocess
import pandas as pd
from scipy.optimize import basinhopping
import copy


class pseudo:
    def __init__(self,file):
        self.file = file
        dr = open(file,'r')
        self.dr = dr

        inv = dict()
        inv_cross = dict()
        uni = dict()
        uni_join = dict()
        are = list()

#        rows_to_skip = 3
        frontmatter = True
        self.tx_section = False
        cnt = 0
        firstentry = True
        find_begin_end = False
        areas = False
        save_before_areas = False
        end = False
        dontrecord = False


        name = ''
        phases_present = list()
        phases_absent = list()
        begin_end = list()
        P = list()
        T = list()
        edges = list()


        for l in dr:

            if frontmatter == True:
                if l[0] == 'x':
                    self.tx_section = True
                    self.tx_section_P = float(l.split()[2])
                if l[0] == 'i' or l[0] == 'u':
                    frontmatter = False

            if frontmatter == False:

                # Use * to locate place in file
                if l[0] == '*' and areas == True:
                    end = True
                if l[0] == '*' and areas == False:
                    areas = True
                    save_before_areas = True

                # Skip rows beginning with space or %
                if l[0] != ' ' and l[0] != '%' and l != '\n':
                    # Skip rows at start of file that describe layout
                    #if cnt >= rows_to_skip:

                    if areas == True and end == False and l[0] != '*':
                        split = l.split()[1:]
                        found_end_bounds = False
                        bounds = list()
                        phases = list()
                        for i in split:
                            if found_end_bounds == False and i !='%':
                                bounds.append(i)
                            if found_end_bounds == True:
                                phases.append(i)
                            if i == '%':
                                found_end_bounds = True
                        are.append([bounds,phases])


                    if areas == False or save_before_areas == True:
                        save_before_areas = False
                        # RECORD EXTRACTED DATA
                        if ((l[0] == 'i' or l[0] == 'u' or l[0] == '*'
                        or l.find('crossover') > -1 or l.find('connect') > -1)):

                            if l.find('crossover') > -1:
                                split = l.split()
                                inv_cross[name] = split[1:] + [phases_present,phases_absent]
                                dontrecord = True
                            if l.find('connect') > -1:
                                split = l.split()
                                uni_join[name] = split[:2] + [phases_present,phases_absent,[split[0],split[1]]]
                                dontrecord = True

                            if dontrecord == False and find_begin_end == False:
                                if firstentry == True:
                                    firstentry = False
                                else:
                                    if name[0] == 'i' and name not in inv_cross.keys():
                                        T = [float(i) for i in T]
                                        P = [float(i) for i in P]
                                        x = [float(i) for i in x]
                                        if self.tx_section == True:
                                            inv[name] = [phases_present,phases_absent,x,T,P]
                                        else:
                                            inv[name] = [phases_present,phases_absent,P[0],T[0]]
                                    if name[0] == 'u' and name not in uni_join.keys():
                                        T = [float(i) for i in T]
                                        P = [float(i) for i in P]
                                        x = [float(i) for i in x]
                                        if self.tx_section == True:
                                            uni[name] = [phases_present,phases_absent,x,T,begin_end,P]
                                        else:
                                            uni[name] = [phases_present,phases_absent,P,T,begin_end]
                            dontrecord = False

                        if find_begin_end == True:
                            split = l.split()
                            begin_end = split[:2]
                            find_begin_end = False

                        # FIND START OF INVARIENT OR UNIVARIENT
                        elif (l[0] == 'i' or l[0] == 'u') and find_begin_end == False:
                            split = l.split()
                            name = split[0]

                            # Extract phases absent and present
                            phases_present = list()
                            phases_absent = list()
                            absentcheck = False
                            for c in split[1:]:
                                if c != '-' and absentcheck == False:
                                    phases_present.append(c)
                                if c == '-':
                                    absentcheck = True
                                if c != '-' and absentcheck == True:
                                    phases_absent.append(c)

                            # Set up P and T lists
                            P = list()
                            T = list()
                            x = list()
                            if name[0] == 'u':
                                find_begin_end = True


                        # RECORD P AND T
                        elif l[0] != '*':
                            split = l.split()
                            if self.tx_section == True:
                                x.append(split[0])
                                P.append(split[1])
                                T.append(split[2])
                            else:
                                P.append(split[0])
                                T.append(split[1])



                    if end == True:
                        split = l.split()
                        if split[0] == 'window':
                            edges = split[1:]

#                    else:
#                        cnt = cnt + 1

        for ptname in inv_cross.keys():
            pt = inv_cross[ptname]
            interp1 = interp1d(uni[pt[0]][2],uni[pt[0]][3])
            interp2 = interp1d(uni[pt[1]][2],uni[pt[1]][3])
            minx = np.max(np.array([np.min(np.array(uni[pt[0]][2])),np.min(np.array(uni[pt[0]][2]))]))
            maxx = np.min(np.array([np.max(np.array(uni[pt[0]][2])),np.max(np.array(uni[pt[0]][2]))]))

            def tosolve(x,interp1=interp1,interp2=interp2,minx=minx,maxx=maxx):
                if x >= minx and x <= maxx:
                    return interp1(x) - interp2(x)
                if x < minx:
                    return (minx-x)**2*(interp1(minx) - interp2(minx))
                if x > maxx:
                    return (maxx-x)**2*(interp1(maxx) - interp2(maxx))

            x0 = newton(tosolve,x0=(maxx-minx)/2+minx)

            inv[ptname] = [pt[2],pt[3],x0,float(interp1(x0))]

        for lnname in uni_join.keys():
            line = uni_join[lnname]
            P = [inv[line[0]][2],inv[line[1]][2]]
            T = [inv[line[0]][3],inv[line[1]][3]]
            uni[lnname] = [line[2],line[3],P,T,line[4][:2]]


        self.inv = inv
        self.inv_cross = inv_cross
        self.uni = uni
        self.uni_join = uni_join
        self.are = are
        if self.tx_section == True:
            self.edges = {'lhs':float(edges[2]),
                          'rhs':float(edges[3]),
                          'bot':float(edges[0]),
                          'top':float(edges[1])}
            self.D_edges = Delaunay([[edges[2],edges[0]],[edges[3],edges[0]],
                                      [edges[2],edges[1]],[edges[3],edges[1]]])
        else:
            self.edges = {'lhs':float(edges[0]),
                          'rhs':float(edges[1]),
                          'bot':float(edges[2]),
                          'top':float(edges[3])}
            self.D_edges = Delaunay([[edges[0],edges[2]],[edges[1],edges[2]],
                                      [edges[0],edges[3]],[edges[1],edges[3]]])

        self.tidylines()
        self.createareas()


    # Plots all raw data crudely
    def crudeplot(self):
        f,a = plt.subplots()
        for line in self.uni.values():
            a.plot(line[3],line[2])
        for pt in self.inv.values():
            a.scatter(pt[3],pt[2],marker='o')
        plt.show()

    # This function makes sure the lines start and end at the points specified
    # in the drawpd file, or start/end at the boundaries of the phase diagram
    # area. Saves output as a new objects- tidyuni.
    def tidylines(self):

        self.tidyuni = dict()

        if self.tx_section == True:
            for invname in self.inv.keys():
                interp_T = interp1d(self.inv[invname][4],self.inv[invname][3])
                interp_x = interp1d(self.inv[invname][4],self.inv[invname][2])

                self.inv[invname][4] = self.tx_section_P
                self.inv[invname][3] = float(interp_T(self.tx_section_P))
                self.inv[invname][2] = float(interp_x(self.tx_section_P))

        for linename in self.uni.keys():
            line = self.uni[linename]
            [begin,end] = line[4]
            T = line[3]
            P = line[2]
            numberpoints = len(T)
            begin_inv_out_of_plot = False
            end_inv_out_of_plot = False
            # Check if line is in the the diagram
            proceed = False
            if any(self.D_edges.find_simplex(np.array([T,P]).T) >= 0):
                proceed = True
            elif begin[0] == 'i':
                if self.D_edges.find_simplex([self.inv[begin][3],self.inv[begin][2]]) >= 0:
                    proceed = True
            elif end[0] == 'i':
                if self.D_edges.find_simplex([self.inv[end][3],self.inv[end][2]]) >= 0:
                    proceed = True
            if proceed == True:
                # Create a list of step changes in T and P
                # Used to identify which variable should be iterated through,
                # i.e. was it calculate for T at P or P at T?
                # Also checks if dT can be used for doing edge interp
                dP = list()
                dT = list()
                for i in range(len(T)-1):
                    dP.append(P[i+1]-P[i])
                    dT.append(T[i+1]-T[i])

                if all(np.array(dP) > 0):
                    useP = True
                else:
                    useP = False
                if all(np.array(dT) > 0):
                    useT = True
                else:
                    useT = False

                if begin[0] == 'i':
                    if self.D_edges.find_simplex([self.inv[begin][3],self.inv[begin][2]]) == -1:
                        begin_inv_out_of_plot = True
                    else:
                        if useP == True:
                            beginval = self.inv[begin][2]
                        else:
                            beginval = self.inv[begin][3]


                if begin == 'begin' or begin_inv_out_of_plot == True:
                    # Find edge it intersects
                    FoundEdge = False
                    i = 0
                    while FoundEdge == False and i<numberpoints:
                        if T[i] <= self.edges['lhs'] and P[i] < self.edges['top'] and P[i] > self.edges['bot']:
                            side = 'lhs'
                            FoundEdge = True
                        elif P[i] >= self.edges['top'] and T[i] < self.edges['rhs'] and T[i] > self.edges['lhs']:
                            side = 'top'
                            FoundEdge = True
                        elif T[i] >= self.edges['rhs'] and P[i] < self.edges['top'] and P[i] > self.edges['bot']:
                            side = 'rhs'
                            FoundEdge = True
                        elif P[i] <= self.edges['top'] and T[i] < self.edges['rhs'] and T[i] > self.edges['lhs']:
                            side = 'bot'
                            FoundEdge = True
                        else:
                            i = i+1


                    if side == 'bot' or side == 'top':
                        if useP == True:
                            beginval = self.edges[side]
                        else:
                            if P[0] == self.edges[side]:
                                beginvalue = T[0]
                            else:
                                if useT == True:
                                    interp = interp1d(P,T)
                                    beginval = interp(float(self.edges[side]))
                                if useT == False:
                                    Tcut = list()
                                    Pcut = list()
                                    if dT[0] > 0:
                                        for i in range(len(T)):
                                            if dT[i] > 0:
                                                Tcut.append(T[i])
                                                Pcut.append(P[i])
                                    if dT[0] < 0:
                                        for i in range(len(T)):
                                            if dT[i] < 0:
                                                Tcut.append(T[i])
                                                Pcut.append(P[i])
                                    interp = interp1d(Pcut,Tcut)
                                    beginval = interp(float(self.edges[side]))


                    if side == 'lhs' or side == 'rhs':
                        # First see if the first point coincides with the edge
                        if useP == True:
                            if T[0] == self.edges[side]:
                                beginval = P[0]
                            else:
                                if useT == True:
                                    interp = interp1d(T,P)
                                    beginval = interp(float(self.edges[side]))
                                if useT == False:
                                    Tcut = list()
                                    Pcut = list()
                                    if dT[0] > 0:
                                        for i in range(len(T)):
                                            if dT[i] > 0:
                                                Tcut.append(T[i])
                                                Pcut.append(P[i])
                                    if dT[0] < 0:
                                        for i in range(len(T)):
                                            if dT[i] < 0:
                                                Tcut.append(T[i])
                                                Pcut.append(P[i])
                                    interp = interp1d(Tcut,Pcut)
                                    beginval = interp(float(self.edges[side]))
                        else:
                            beginval = self.edges[side]

                if end[0] == 'i':
                    if self.D_edges.find_simplex([self.inv[end][3],self.inv[end][2]]) == -1:
                        end_inv_out_of_plot = True
                    else:
                        if useP == True:
                            endval = self.inv[end][2]
                        else:
                            endval = self.inv[end][3]


                if end == 'end' or end_inv_out_of_plot == True:
                    # Find edge it intersects
                    FoundEdge = False
                    i = 1
                    while FoundEdge == False and i <= numberpoints:
                        if T[-i] <= self.edges['lhs'] and P[-i] < self.edges['top'] and P[-i] > self.edges['bot']:
                            side = 'lhs'
                            FoundEdge = True
                        elif P[-i] >= self.edges['top'] and T[-i] < self.edges['rhs'] and T[-i] > self.edges['lhs']:
                            side = 'top'
                            FoundEdge = True
                        elif T[-i] >= self.edges['rhs'] and P[-i] < self.edges['top'] and P[-i] > self.edges['bot']:
                            side = 'rhs'
                            FoundEdge = True
                        elif P[-i] <= self.edges['top'] and T[-i] < self.edges['rhs'] and T[-i] > self.edges['lhs']:
                            side = 'bot'
                            FoundEdge = True
                        else:
                            i = i+1

                    if side == 'bot' or side == 'top':
                        if useP == True:
                            endval = self.edges[side]
                        else:
                            if P[0] == self.edges[side]:
                                beginvalue = T[0]
                            else:
                                if useT == True:
                                    interp = interp1d(P,T)
                                    beginval = interp(float(self.edges[side]))
                                if useT == False:
                                    Tcut = list()
                                    Pcut = list()
                                    if dT[0] > 0:
                                        for i in range(len(T)-1):
                                            if dT[i] > 0:
                                                Tcut.append(T[i+1])
                                                Pcut.append(P[i+1])
                                    if dT[0] < 0:
                                        for i in range(len(T)-1):
                                            if dT[i] < 0:
                                                Tcut.append(T[i+1])
                                                Pcut.append(P[i+1])
                                    interp = interp1d(Pcut,Tcut)
                                    beginval = interp(float(self.edges[side]))

                    if side == 'lhs' or side == 'rhs':
                        if useP == True:
                            if T[-1] == self.edges[side]:
                                endval = P[-1]
                            else:
                                if useT == True:
                                    interp = interp1d(T,P)
                                    endval = interp(float(self.edges[side]))
                                if useT == False:
                                    Tcut = list()
                                    Pcut = list()
                                    if dT[-1] > 0:
                                        for i in range(len(T)-1):
                                            if dT[i] > 0:
                                                Tcut.append(T[i])
                                                Pcut.append(P[i])
                                    if dT[-1] < 0:
                                        for i in range(len(T)-1):
                                            if dT[i] < 0:
                                                Tcut.append(T[i])
                                                Pcut.append(P[i])
                                    interp = interp1d(Tcut,Pcut)
                                    endval = interp(float(self.edges[side]))
                        else:
                            beginval = self.edges[side]


                # Construct new P, T lists
                newT = list()
                newP = list()


                interpP = interp1d(P,T,bounds_error=False,fill_value='extrapolate')
                interpT = interp1d(T,P,bounds_error=False,fill_value='extrapolate')

                if begin[0] == 'i' and begin_inv_out_of_plot == False:
                    newP.append(self.inv[begin][2])
                    newT.append(self.inv[begin][3])
                else:
                    if useP == True:
                        newP.append(beginval)
                        newT.append(interpP(beginval))
                    else:
                        newT.append(beginval)
                        newP.append(interpT(beginval))

                if linename not in self.uni_join.keys():
                    for i in range(len(T)):
                        if useP == True:
                            if ((T[i] > interpP(beginval) and T[i] < interpP(endval)) or (P[i] > beginval and P[i] < endval))\
                                    or ((T[i] < interpP(beginval) and T[i] > interpP(endval)) or (P[i] < beginval and P[i] > endval)):
                                newP.append(P[i])
                                newT.append(T[i])
                        else:
                            if ((T[i] > beginval and T[i] < endval) or (P[i] > interpT(beginval) and P[i] < interpT(endval)))\
                                    or ((T[i] < beginval and T[i] > endval) or (P[i] < interpT(beginval) and P[i] > interpT(endval))):
                                newP.append(P[i])
                                newT.append(T[i])


                if end[0] == 'i' and end_inv_out_of_plot == False:
                    newP.append(self.inv[end][2])
                    newT.append(self.inv[end][3])
                else:
#                    if endval > np.min(P) and endval < np.max(P):
                    if useP == True:
                        newP.append(endval)
                        newT.append(interpP(endval))
                    else:
                        newP.append(interpT(endval))
                        newT.append(endval)

                self.tidyuni[linename] = [line[0],line[1],newP,newT]

    def createareas(self):
        self.tidyarea = list()
        self.polygons = list()
        # Calculate minimum and maximum phases in areas on phase diagram
        numphases = list()
        for area in self.are:
            numphases.append(len(area[1]))
        numphases = np.array(numphases)
        self.minphases = np.min(numphases)
        self.maxphases = np.max(numphases)


        # First work out which areas occupy the corners of the phase diagram.

        # Find the lines which intersect the sides
        rhs_lines = list()
        lhs_lines = list()
        top_lines = list()
        bot_lines = list()
        rhs_lines_n = list()
        lhs_lines_n = list()
        top_lines_n = list()
        bot_lines_n = list()

        for linename in self.tidyuni.keys():
            line = self.tidyuni[linename]
            if line[3][0] == self.edges['lhs']:
                lhs_lines.append(line[2][0])
                lhs_lines_n.append(linename)
            if line[3][-1] == self.edges['lhs']:
                lhs_lines.append(line[2][-1])
                lhs_lines_n.append(linename)

            if line[3][0] == self.edges['rhs']:
                rhs_lines.append(line[2][0])
                rhs_lines_n.append(linename)
            if line[3][-1] == self.edges['rhs']:
                rhs_lines.append(line[2][-1])
                rhs_lines_n.append(linename)

            if line[2][0] == self.edges['top']:
                top_lines.append(line[3][0])
                top_lines_n.append(linename)
            if line[2][-1] == self.edges['top']:
                top_lines.append(line[3][-1])
                top_lines_n.append(linename)

            if line[2][0] == self.edges['bot']:
                bot_lines.append(line[3][0])
                bot_lines_n.append(linename)
            if line[2][-1] == self.edges['bot']:
                bot_lines.append(line[3][-1])
                bot_lines_n.append(linename)


        # Find the lines closest to the corners, or does the line span the diagram?
        ul = [False,False]
        ur = [False,False]
        ll = [False,False]
        lr = [False,False]
        if len(lhs_lines) > 0:
            ul[0] = lhs_lines_n[np.argmax(lhs_lines)]
            ll[1] = lhs_lines_n[np.argmin(lhs_lines)]
        if len(top_lines) > 0:
            ul[1] = top_lines_n[np.argmin(top_lines)]
            ur[0] = top_lines_n[np.argmax(top_lines)]
        if len(rhs_lines) > 0:
            ur[1] = rhs_lines_n[np.argmax(rhs_lines)]
            lr[0] = rhs_lines_n[np.argmin(rhs_lines)]
        if len(bot_lines) > 0:
            lr[1] = bot_lines_n[np.argmax(bot_lines)]
            ll[0] = bot_lines_n[np.argmin(bot_lines)]


        # Process each area
        for area in self.are:
            areaP = list()
            areaT = list()
            # Check if the area contains a corner

            # Check if it is a corner rather than a side
            cornerline = False
            foundcorner = False
            if ul[0] != False and ul[1] != False:
                if ul[0] in area[0] and ul[1] in area[0]:
                    #if ul[0] != ul[1]:
                    #    foundcorner = True
                    if (ul[0] == ul[1] and len(area[0]) == 1) or ul[0] != ul[1]:
                        foundcorner = True
                        cornerline = ul[0]
                        if self.tidyuni[ul[0]][3][0] == self.edges['lhs']:
                            cornerlineT = [self.edges['lhs']] + self.tidyuni[ul[0]][3]
                            cornerlineP = [self.edges['top']] + self.tidyuni[ul[0]][2]
                        else:
                            cornerlineT = self.tidyuni[ul[0]][3] + [self.edges['lhs']]
                            cornerlineP = self.tidyuni[ul[0]][2] + [self.edges['top']]
            if ur[0] != False and ur[1] != False and foundcorner != True:
                if ur[0] in area[0] and ur[1] in area[0]:
                    #if ul[0] != ul[1]:
                    #    foundcorner = True
                    if (ur[0] == ur[1] and len(area[0]) == 1) or ur[0] != ur[1]:
                        foundcorner = True
                        cornerline = ur[0]
                        if self.tidyuni[ur[0]][2][0] == self.edges['top']:
                            cornerlineT = [self.edges['rhs']] + self.tidyuni[ur[0]][3]
                            cornerlineP = [self.edges['top']] + self.tidyuni[ur[0]][2]
                        else:
                            cornerlineT = self.tidyuni[ur[0]][3] + [self.edges['rhs']]
                            cornerlineP = self.tidyuni[ur[0]][2] + [self.edges['top']]
            if lr[0] != False and lr[1] != False and foundcorner != True:
                if lr[0] in area[0] and lr[1] in area[0]:
                    #if ul[0] != ul[1]:
                    #    foundcorner = True
                    if (lr[0] == lr[1] and len(area[0]) == 1) or lr[0] != lr[1]:
                        foundcorner = True
                        cornerline = lr[0]
                        if self.tidyuni[lr[0]][3][0] == self.edges['rhs']:
                            cornerlineT = [self.edges['rhs']] + self.tidyuni[lr[0]][3]
                            cornerlineP = [self.edges['bot']] + self.tidyuni[lr[0]][2]
                        else:
                            cornerlineT = self.tidyuni[lr[0]][3] + [self.edges['rhs']]
                            cornerlineP = self.tidyuni[lr[0]][2] + [self.edges['bot']]
            if ll[0] != False and ll[1] != False and foundcorner != True:
                if ll[0] in area[0] and ll[1] in area[0]:
                    #if ul[0] != ul[1]:
                    #foundcorner = True
                    if (ll[0] == ll[1] and len(area[0]) == 1) or ll[0] != ll[1]:
                        foundcorner = True
                        cornerline = ll[0]
                        if self.tidyuni[ll[0]][2][0] == self.edges['bot']:
                            cornerlineT = [self.edges['lhs']] + self.tidyuni[ll[0]][3]
                            cornerlineP = [self.edges['bot']] + self.tidyuni[ll[0]][2]
                        else:
                            cornerlineT = self.tidyuni[ll[0]][3] + [self.edges['lhs']]
                            cornerlineP = self.tidyuni[ll[0]][2] + [self.edges['bot']]


            # Check if there's a side
            cornerlines = [None,None]
            if ul[0] == False and ll[1] == False:
                if (ul[1] in area[0] and ll[0] in area[0]):
                    cornerlines = [ul[1],ll[0]]
                    if self.tidyuni[ul[1]][2][0] == self.edges['top']:
                        cornerlineT = self.tidyuni[ul[1]][3][::-1] + [self.edges['lhs']]
                        cornerlineP = self.tidyuni[ul[1]][2][::-1] + [self.edges['top']]
                    else:
                        cornerlineT = self.tidyuni[ul[1]][3] + [self.edges['lhs']]
                        cornerlineP = self.tidyuni[ul[1]][2] + [self.edges['top']]

                    if self.tidyuni[ll[0]][2][0] == self.edges['bot']:
                        cornerlineT += [self.edges['lhs']] + self.tidyuni[ll[0]][3]
                        cornerlineP += [self.edges['bot']] + self.tidyuni[ll[0]][2]
                    else:
                        cornerlineT += [self.edges['lhs']] + self.tidyuni[ll[0]][3][::-1]
                        cornerlineP += [self.edges['bot']] + self.tidyuni[ll[0]][2][::-1]


            if ul[1] == False and ur[0] == False:
                if (ul[0] in area[0] and ur[1] in area[0]):
                    cornerlines = [ul[0],ur[1]]
                    if self.tidyuni[ul[0]][3][0] == self.edges['lhs']:
                        cornerlineT = self.tidyuni[ul[0]][3][::-1] + [self.edges['lhs']]
                        cornerlineP = self.tidyuni[ul[0]][2][::-1] + [self.edges['top']]
                    else:
                        cornerlineT = self.tidyuni[ul[0]][3] + [self.edges['lhs']]
                        cornerlineP = self.tidyuni[ul[0]][2] + [self.edges['top']]
                    if self.tidyuni[ur[1]][3][0] == self.edges['rhs']:
                        cornerlineT += [self.edges['rhs']] + self.tidyuni[ur[1]][3]
                        cornerlineP += [self.edges['top']] + self.tidyuni[ur[1]][2]
                    else:
                        cornerlineT += [self.edges['rhs']] + self.tidyuni[ur[1]][3][::-1]
                        cornerlineP += [self.edges['top']] + self.tidyuni[ur[1]][2][::-1]


            if ur[1] == False and lr[0] == False:
                if (ur[0] in area[0] and lr[1] in area[0]):
                    cornerlines = [ur[0],lr[1]]
                    if self.tidyuni[ur[0]][2][0] == self.edges['top']:
                        cornerlineT = self.tidyuni[ur[0]][3][::-1] + [self.edges['rhs']]
                        cornerlineP = self.tidyuni[ur[0]][2][::-1] + [self.edges['top']]
                    else:
                        cornerlineT = self.tidyuni[ur[0]][3] + [self.edges['rhs']]
                        cornerlineP = self.tidyuni[ur[0]][2] + [self.edges['top']]
                    if self.tidyuni[lr[1]][2][0] == self.edges['bot']:
                        cornerlineT += [self.edges['rhs']] + self.tidyuni[lr[1]][3]
                        cornerlineP += [self.edges['bot']] + self.tidyuni[lr[1]][2]
                    else:
                        cornerlineT += [self.edges['rhs']] + self.tidyuni[lr[1]][3][::-1]
                        cornerlineP += [self.edges['bot']] + self.tidyuni[lr[1]][2][::-1]


            if lr[1] == False and ll[0] == False:
                if (lr[0] in area[0] and ll[1] in area[0]):
                    cornerlines = [lr[0],ll[1]]
                    if self.tidyuni[lr[0]][3][0] == self.edges['rhs']:
                        cornerlineT = self.tidyuni[lr[0]][3][::-1] + [self.edges['rhs']]
                        cornerlineP = self.tidyuni[lr[0]][2][::-1] + [self.edges['bot']]
                    else:
                        cornerlineT = self.tidyuni[lr[0]][3] + [self.edges['rhs']]
                        cornerlineP = self.tidyuni[lr[0]][2] + [self.edges['bot']]
                    if self.tidyuni[ll[1]][2][0] == self.edges['lhs']:
                        cornerlineT += [self.edges['lhs']] + self.tidyuni[ll[1]][3]
                        cornerlineT += [self.edges['bot']] + self.tidyuni[ll[1]][2]
                    else:
                        cornerlineT += [self.edges['lhs']] + self.tidyuni[ll[1]][3][::-1]
                        cornerlineT += [self.edges['bot']] + self.tidyuni[ll[1]][2][::-1]



            firstline = True

            for linename in area[0]:
                if cornerlines != [None,None]:
                    if linename == cornerlines[0]:
                        lineP = cornerlineP
                        lineT = cornerlineT

                    #elif linename != cornerlines[1]:
                    elif linename != cornerlines[1]:
                        lineP = self.tidyuni[linename][2]
                        lineT = self.tidyuni[linename][3]

                else:
                    if linename == cornerline:
                        lineP = cornerlineP
                        lineT = cornerlineT
                    else:
                        lineP = self.tidyuni[linename][2]
                        lineT = self.tidyuni[linename][3]

                if firstline == False:
                    if areaP[-1] == lineP[0] or areaT[-1] == lineT[0]:
                        areaP = areaP + lineP
                        areaT = areaT + lineT
                    elif areaP[0] == lineP[0] or areaT[0] == lineT[0]:
                        areaP = areaP[::-1] + lineP
                        areaT = areaT[::-1] + lineT
                    elif areaP[-1] == lineP[-1] or areaT[-1] == lineT[-1]:
                        areaP = areaP + lineP[::-1]
                        areaT = areaT + lineT[::-1]
                    elif areaP[0] == lineP[-1] or areaT[0] == lineT[-1]:
                        areaP = areaP[::-1] + lineP[::-1]
                        areaT = areaT[::-1] + lineT[::-1]

                elif linename != cornerlines[1] or cornerlines[0] == cornerlines[1]:
                    areaP = lineP
                    areaT = lineT
                    firstline = False


            self.tidyarea.append([areaP,areaT,area[1]])
            self.polygons.append([Polygon(tuple(map(tuple,np.array([areaP,areaT]).T))),area[1]])





    def imitate_drawpd(self,color=(0.2,0.0,1.0)):
        f,a = plt.subplots(figsize=[5,5])

        if self.tx_section == False:
            for area in self.tidyarea:
                c = len(area[2])
                c = (c-self.minphases)/(self.maxphases-self.minphases)
                c = (color[0] + (1-color[0])*c, color[1]+ (1-color[1])*c, color[2]+ (1-color[2])*c)

                a.fill(area[1],area[0],color=c)

            for line in self.tidyuni.values():
                a.plot(line[3],line[2],c='k',lw=0.8)

            a.set_xlim(self.edges['lhs'],self.edges['rhs'])
            a.set_ylim(self.edges['bot'],self.edges['top'])

            a.set_xlabel('Temperature ($^\circ$C)')
            a.set_ylabel('Pressure (kbar)')
        else:
            for area in self.tidyarea:
                c = len(area[2])
                c = (c-self.minphases)/(self.maxphases-self.minphases)
                c = (color[0] + (1-color[0])*c, color[1]+ (1-color[1])*c, color[2]+ (1-color[2])*c)

                a.fill(area[0],area[1],color=c)

            for line in self.tidyuni.values():
                a.plot(line[2],line[3],c='k',lw=0.8)

            a.set_ylim(self.edges['lhs'],self.edges['rhs'])
            a.set_xlim(self.edges['bot'],self.edges['top'])

            a.set_ylabel('Temperature ($^\circ$C)')
            a.set_xlabel('x')


        return f,a

    def makegrid(self,gridT,gridP):
        gridT, gridP = np.meshgrid(gridT,gridP)

        phases = np.full([np.shape(gridT)[0],np.shape(gridT)[1],self.maxphases],'   ')

        for x in range(np.shape(gridT)[0]):
            for y in range(np.shape(gridT)[1]):
                foundarea = False
                i=0
                while foundarea == False and i < len(self.tidyarea):
                    if self.polygons[i][0].contains(Point((gridP[x,y],gridT[x,y]))) == True:
                        foundarea = True
                        phases[x,y,:len(self.tidyarea[i][-1])] = self.tidyarea[i][-1]
                    else:
                        i=i+1

        return gridT, gridP, phases


class tcgridrun:
    def __init__(self,gridT,gridP,phases,tx_section=False,tx_section_P=50.0,onlygaps=False,results=None):

        if tx_section == False:
            self.gridT = gridT
            self.gridP = gridP
            self.phases = phases

            self.commands = b''

            for x in range(np.shape(gridT)[0]):
                for y in range(np.shape(gridT)[1]):
                    run_tc = True
                    if onlygaps == True:
                        # Need to check if both exist in same record!!!!
                        if results[(results.temperature==gridT[x,y])&(results.pressure==gridP[x,y])].shape[0] == 1:
                            run_tc = False

                    if run_tc == True:
                        for p in self.phases[x,y,:]:
                            if p != '   ':
                                self.commands = self.commands + p.encode() + b' '
                        self.commands = self.commands + b'\n \n \n'

                        self.commands = self.commands + str(gridP[x,y]).encode() + b' ' + str(gridP[x,y]).encode() + b' \n '
                        self.commands = self.commands + str(gridT[x,y]).encode() + b' ' + str(gridT[x,y]).encode() + b' \n '

                        if x == np.shape(gridT)[0]-1 and y == np.shape(gridT)[1]-1:
                            self.commands = self.commands + b' n \n \n '
                        else:
                            self.commands = self.commands + b' y \n '

        else:
            self.gridT = gridT
            self.gridP = gridP
            self.phases = phases
            self.tx_section_P = tx_section_P

            self.commands = b''

            for x in range(np.shape(gridT)[0]):
                for y in range(np.shape(gridT)[1]):
                    run_tc = True
                    if onlygaps == True:
                        # Need to check if both exist in same record!!!!
                        if results[(results.temperature==gridT[x,y])&(results.x==np.round(gridP[x,y],5))].shape[0] > 0:
                            run_tc = False

                    if run_tc == True:
                        for p in self.phases[x,y,:]:
                            if p != '   ':
                                self.commands = self.commands + p.encode() + b' '
                        self.commands = self.commands + b'\n \n \n'

                        self.commands = self.commands + str(self.tx_section_P).encode() + b' \n '
                        self.commands = self.commands + str(gridT[x,y]).encode()  + b' \n '
                        self.commands = self.commands + str(np.round(gridP[x,y],5)).encode()  + b' \n '


                        if x == np.shape(gridT)[0]-1 and y == np.shape(gridT)[1]-1:
                            self.commands = self.commands + b' n \n \n '
                        else:
                            self.commands = self.commands + b' y \n '

    def run(self):

        proc = subprocess.Popen("./tc340i",shell=True,
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE)
        proc.communicate(self.commands)
        proc.terminate()

def import_logfile(logfile,save=True,prefix='tc',tx_section=False):
    lf = open(logfile,'r',encoding='mac-roman')

    if tx_section == True:
        results = pd.DataFrame(columns=['pressure','temperature','x'])

        new_entry = False
        first_entry = True
        read_vals = False
        read_mode = False
        read_bc = False
        read_liq_comp = False
        check_next_row = False
        read_overflow_vals = False
        phases = []
        for l in lf:
            split = l.split()

            if len(split) > 0 and l != '\n':

                replace_split = False
                for i,segment in zip(range(len(split)),split):
                    if len(segment.split('.')) > 2:
                        first_stop = segment.find('.')
                        second_stop = segment[first_stop+1:].find('.')
                        third_stop = segment[first_stop+second_stop+2:].find('.')
                        new_segment1 = segment[:first_stop+second_stop]
                        if third_stop != -1:
                            new_segment2 = segment[first_stop+second_stop:first_stop+second_stop+third_stop]
                            new_segment3 = segment[first_stop+second_stop+third_stop:]

                        else:
                            new_segment2 = segment[first_stop+second_stop:]
                        problem_segment = i
                        replace_split = True
                if replace_split == True:
                    new_split = split[:problem_segment]
                    if third_stop != -1:
                        new_split += [new_segment1,new_segment2,new_segment3]
                    else:
                        new_split += [new_segment1,new_segment2]
                    new_split += split[problem_segment+1:]
                    split = new_split


                if read_overflow_vals == True:
                    for i in range(len(split)):
                        row[col_names[i]] = float(split[i])
                    read_overflow_vals = False
                if check_next_row == True:
                    check_next_row = False
                    if split[0] != 'SiO2' and split[0] != 'mode':
                        col_names = split
                        read_overflow_vals = True
                        check_next_row = False
                if l[0] == '^':
                    new_entry = True
                    if first_entry == False and solution_found == True:
                        if 'liq' in list(row.index):
                            row['X'] = MassFracFromMols(bc,lc,row['liq'])
                        else:
                            row['X'] = 0.0
                        results = results.append(row, ignore_index = True)
                    else:
                        first_entry = False
                    row = pd.Series()
                    bc = pd.Series()
                    lc = pd.Series()
                    solution_found = False

                elif split[0] == 'start':
                    row['x'] = float(split[-1])

                elif read_vals == True:
                    check_next_row = True
                    read_vals = False
                    vals = split
                    for i in range(len(vals)):
                        row[col_names[i]] = float(vals[i])
                elif split[0] == 'P(kbar)':
                    col_names = ['pressure','temperature'] + split[2:]
                    read_vals = True
                elif split[0] == 'phases':
                    for p in split[2:]:
                        if p not in phases:
                            phases += [p]
                elif read_mode == True:
                    read_mode = False
                    for i in range(len(col_names)):
                        row[col_names[i]] = float(split[i])
                elif split[0] == 'mode':
                    read_mode = True
                    col_names = split[1:]
                    solution_found = True
                    check_next_row = False
                elif read_bc == True:
                    if split[0] == 'SiO2':
                        col_names = split
                        majors = split
                    else:
                        for i in range(len(col_names)):
                            row['bc_'+col_names[i]] = float(split[i])
                            bc[col_names[i]] = float(split[i])
                            read_bc = False
                elif split[0] == 'composition':
                    read_bc = True
                elif read_liq_comp == True:
                    for i in range(len(majors)):
                        row['liq_'+majors[i]] = float(split[i])
                        read_liq_comp = False
                        lc[majors[i]] = float(split[i])
                elif split[0] == 'SiO2' and read_bc == False:
                    read_liq_comp = True
        if solution_found == True:
            if 'liq' in list(row.index):
                row['X'] = MassFracFromMols(bc,lc,row['liq'])
            else:
                row['X'] = 0.0
            results = results.append(row, ignore_index = True)


        if save == True:
            results.to_csv(prefix+'_tcgrid.csv')
            np.save(prefix+'_phases',np.array(phases))

        return tcgridresults(results,phases,None,None,tx_section=tx_section)

    else:
        print('Version of code for pseudosections is not written.')


def import_icfile(icfile,save=True,prefix='tc'):
    ic = open(icfile,'r',encoding='mac-roman')

    results = pd.DataFrame(columns=['pressure','temperature'])


    results = pd.DataFrame(columns=['pressure','temperature'])
    foundT = False
    foundSpace = 0
    foundbulkcomp = False
    readbulkcomp = False
    readcomp = False
    readMode = False
    extract = {}
    majormols = {}
    phase = 'none'
    phases = []
    for l in ic:
        split = l.split()
        if foundT == True:
            if l == '\n':
                foundSpace = foundSpace + 1

        if len(split)>0:
            if readMode == True:
                if split[0] == 'G':
                    readMode = False
                else:
                    extract[split[0]] = float(split[1])
                    if (split[0] in phases) == False:
                        phases = phases + [split[0]]

            if split[0] == 'mode':
                readMode = True
                readcomp = False

            if readcomp == True:
                if split[0] == 'liq':
                    liqmols = split[1:]
                    tot = 0
                    liqc = {}
                    for i in range(len(majormols)):
                        extract['liq_'+majormols[i]] = np.float(liqmols[i])
                        tot = tot + np.float(liqmols[i])
                    for i in range(len(majormols)):
                        extract['liq_'+majormols[i]] = np.float(liqmols[i])/tot*100
                        liqc[majormols[i]] = np.float(liqmols[i])/tot*100
                elif split[0] != 'bulk':
                    phasemols = split[1:]
                    tot = 0
                    for i in range(len(majormols)):
                        extract[split[0]+'_'+majormols[i]] = np.float(phasemols[i])
                        tot = tot + np.float(phasemols[i])
                    for i in range(len(majormols)):
                        extract[split[0]+'_'+majormols[i]] = np.float(phasemols[i])/tot*100

#                        lookforliqcomp = False

            if readbulkcomp == True:
                bcmols = split[1:]
                readbulkcomp = False
                bc = {}
                for i in range(len(majormols)):
                    bc[majormols[i]] = float(bcmols[i])
                bulkcomp = bc
                bulkcompraw = copy.deepcopy(bc)
                readbulkcomp = False
                foundbulkcomp = True
            if split[0] == 'SiO2':
                if foundbulkcomp == False:
                    readbulkcomp = True
                    majormols = split
                readcomp = True
                foundSpace = 0
                foundT = False
#                        results = results.append(extract,ignore_index=True)
            if split[0] == 'check':
                if 'liq' in extract.keys():
                    extract['X'] = MassFracFromMols(bc,liqc,extract['liq'])
                results = results.append(extract,ignore_index=True)
            if foundSpace >= 2:
                if len(split) == 3:
                    phase = split[0]
                    extract[phase+'_'+split[1]] = split[2]
                elif phase != 'none':
                    extract[phase+'_'+split[0]] = split[1]




            if split[0] == 'P(kbar)':
                extract = {}
                extract['pressure'] = np.round(float(split[1]),2)
            if split[0][0] == 'T':
                extract['temperature'] = np.round(float(split[1]),2)
                foundT = True

    if save == True:
        results.to_csv(prefix+'_tcgrid.csv')
        np.save(prefix+'_phases',np.array(phases))
        pd.Series(bulkcomp).to_csv(prefix+'_bc.csv')
        pd.Series(bulkcompraw).to_csv(prefix+'_bcr.csv')

    return tcgridresults(results,phases,bulkcomp,bulkcompraw)

#        self.results = results

def load_tcgrid(prefix='tc',load_bc=True,tx_section=False):
    results = pd.read_csv(prefix+'_tcgrid.csv')
    phases = np.load(prefix+'_phases.npy')
    if load_bc == True:
        bulkcomp = pd.read_csv(prefix+'_bc.csv',index_col=0,header=None).iloc[:,0]
        bulkcompraw = pd.read_csv(prefix+'_bcr.csv',index_col=0,header=None).iloc[:,0]

        return tcgridresults(results,phases,bulkcomp,bulkcompraw)
    else:
        return tcgridresults(results,phases,None,None,tx_section=tx_section)

class tcgridresults:
    def __init__(self,table,phases,bulkcomp,bulkcompraw,tx_section=False):
        self.results = table
        self.phases = phases
        self.bulkcomp = bulkcomp
        self.bulkcompraw = bulkcompraw
        self.tx_section = tx_section


    def gridresults(self,variable='X'):

        if self.tx_section == False:
            P = np.sort(self.results.pressure.unique())
            T = np.sort(self.results.temperature.unique())

            #X = np.zeros([len(P),len(T)])
            X = np.full([len(P),len(T)],np.nan)

            for x in range(np.shape(X)[0]):
                for y in range(np.shape(X)[1]):
                    if len(self.results[variable][(self.results.pressure==P[x])&(self.results.temperature==T[y])]) > 0:
                        X[x,y] = self.results[variable][(self.results.pressure==P[x])&(self.results.temperature==T[y])]

        else:
            c = np.sort(self.results.x.unique())
            T = np.sort(self.results.temperature.unique())

            X = np.zeros([len(c),len(T)])

            for x in range(np.shape(X)[0]):
                for y in range(np.shape(X)[1]):
                    if len(self.results[variable][(self.results.x==c[x])&(self.results.temperature==T[y])]) > 0:
                        X[x,y] = self.results[variable][(self.results.x==c[x])&(self.results.temperature==T[y])]

        return X

    def CalculateSites(self):
        if 'cpx' in self.results.columns:
            cpx_ends = {'crdi':[],'cess':[],'jd':[],'cats':[],'di':[],'cen':[],'cfm':[],'fs':[]}
            cpx_sites = {'xFeM1':[],'xFe3M1':[],'xFeM2':[],'FeO':[]}
        if 'opx' in self.results.columns:
            opx_ends = {'mess':[],'cren':[],'odi':[],'mgts':[],'fs':[],'fm':[],'en':[]}
            opx_sites = {'xFeM1':[],'xFe3M1':[],'xFeM2':[],'FeO':[]}
        if 'ol' in self.results.columns:
            ol_ends = {'fa':[],'fo':[]}
            ol_sites = {'xFeM':[],'FeO':[]}
        if 'g' in self.results.columns:
            g_ends = {'andr':[],'knor':[],'gr':[],'alm':[],'py':[]}
            g_sites = {'xFeM1':[],'xFe3M2':[],'FeO':[]}

        for i, row in self.results.iterrows():
            if 'cpx' in self.results.columns and row['cpx'] > 0:
                cpx_ends['crdi'].append(row['cr(cpx)'])
                cpx_ends['cess'].append(row['f(cpx)'])
                cpx_ends['jd'].append(row['n(cpx)'])
                cpx_ends['cats'].append(row['y(cpx)']-row['cr(cpx)']-row['f(cpx)'])
                cpx_ends['di'].append(1.0 - row['o(cpx)'] - row['n(cpx)'] - row['y(cpx)'])
                cpx_ends['cen'].append(row['o(cpx)']*(1.0-row['x(cpx)']) + row['Q(cpx)']*
                                        (1-row['y(cpx)']-row['n(cpx)']))
                cpx_ends['cfm'].append((row['n(cpx)']+row['y(cpx)']-1)*(2*row['Q(cpx)']+row['x(cpx)']) +
                                        row['o(cpx)']*row['x(cpx)'])
                cpx_ends['fs'].append((row['Q(cpx)']+row['x(cpx)'])*(1-row['y(cpx)']-row['n(cpx)']))

                cpx_sites['xFeM1'].append(cpx_ends['fs'][-1])
                cpx_sites['xFeM2'].append(cpx_ends['fs'][-1]+cpx_ends['cfm'][-1])
                cpx_sites['xFe3M1'].append(cpx_ends['cess'][-1])
                cpx_sites['FeO'].append(50*(cpx_sites['xFeM1'][-1]+cpx_sites['xFeM2'][-1]+cpx_sites['xFe3M1'][-1]))
            elif 'cpx' in self.results.columns:
                for k in list(cpx_ends.keys()):
                    cpx_ends[k].append(0.0)
                for k in list(cpx_sites.keys()):
                    cpx_sites[k].append(0.0)
            if 'opx' in self.results.columns and row['opx'] > 0:
                opx_ends['mess'].append(row['f(opx)'])
                opx_ends['cren'].append(row['cr(opx)'])
                opx_ends['odi'].append(row['c(opx)'])
                opx_ends['mgts'].append(row['y(opx)']-row['cr(opx)']-row['f(opx)'])
                opx_ends['fs'].append((row['Q(opx)']+row['x(opx)'])*(1.0-row['y(opx)']))
                opx_ends['fm'].append((2*row['Q(opx)']+row['x(opx)'])*(row['y(opx)']-1.0))
                opx_ends['en'].append((row['Q(opx)']+1.0)*(1.0-row['y(opx)'])-row['c(opx)'])

                opx_sites['xFeM1'].append(opx_ends['fs'][-1])
                opx_sites['xFe3M1'].append(opx_ends['mess'][-1])
                opx_sites['xFeM2'].append(opx_ends['fs'][-1]+opx_ends['fm'][-1])
                opx_sites['FeO'].append(50*(opx_sites['xFeM1'][-1]+opx_sites['xFeM2'][-1]+opx_sites['xFe3M1'][-1]))
            elif 'opx' in self.results.columns:
                for k in list(opx_ends.keys()):
                    opx_ends[k].append(0.0)
                for k in list(opx_sites.keys()):
                    opx_sites[k].append(0.0)
            if 'ol' in self.results.columns and row['ol'] > 0:
                ol_ends['fo'].append(1.0-row['x(ol)'])
                ol_ends['fa'].append(row['x(ol)'])

                ol_sites['xFeM'].append(ol_ends['fa'][-1])
                ol_sites['FeO'].append(ol_sites['xFeM'][-1]*100)
            elif 'ol' in self.results.columns:
                for k in list(ol_ends.keys()):
                    ol_ends[k].append(0.0)
                for k in list(ol_sites.keys()):
                    ol_sites[k].append(0.0)
            if 'g' in self.results.columns and row['g'] > 0:
                g_ends['andr'].append(row['f(g)'])
                g_ends['knor'].append(row['cr(g)'])
                g_ends['gr'].append(row['c(g)']-row['f(g)'])
                g_ends['alm'].append(row['x(g)']*(1-row['c(g)']))
                g_ends['py'].append((1-row['x(g)'])*(1-row['c(g)'])-row['cr(g)'])

                g_sites['xFeM1'].append(g_ends['alm'][-1])
                g_sites['xFe3M2'].append(g_ends['andr'][-1])
                g_sites['FeO'].append(100*(0.6*g_sites['xFeM1'][-1]+0.4*g_sites['xFe3M2'][-1]))
            elif 'g' in self.results.columns:
                for k in list(g_ends.keys()):
                    g_ends[k].append(0.0)
                for k in list(g_sites.keys()):
                    g_sites[k].append(0.0)

        if 'cpx' in self.results.columns:
            for k in list(cpx_ends.keys()):
                self.results['cpx_'+k] = cpx_ends[k]
            for k in list(cpx_sites.keys()):
                self.results['cpx_'+k] = cpx_sites[k]
        if 'opx' in self.results.columns:
            for k in list(opx_ends.keys()):
                self.results['opx_'+k] = opx_ends[k]
            for k in list(opx_sites.keys()):
                self.results['opx_'+k] = opx_sites[k]
        if 'ol' in self.results.columns:
            for k in list(ol_ends.keys()):
                self.results['ol_'+k] = ol_ends[k]
            for k in list(ol_sites.keys()):
                self.results['ol_'+k] = ol_sites[k]
        if 'g' in self.results.columns:
            for k in list(g_ends.keys()):
                self.results['g_'+k] = g_ends[k]
            for k in list(g_sites.keys()):
                self.results['g_'+k] = g_sites[k]







def MassFracFromMols(bc,pc,f):
    oxidemass = {'SiO2':60.0855,
                 'TiO2':79.867,
                 'Al2O3':101.964,
                 'FeO':71.845,
                 'Fe2O3':159.69,
                 'MnO':70.938,
                 'MgO':40.305,
                 'CaO':56.078,
                 'Na2O':61.98,
                 'K2O':94.2,
                 'Cr2O3':150}
#    elementmass = {'Si':28.0855,
#                   'Ti':47.867,
#                   'Al':26.982,
#                   'Fe':55.845,
#                   'Mn':54.938,
#                   }

    if 'O' in pc.keys():
        pc['Fe2O3'] = pc['O']
        pc['FeO'] = pc['FeO'] - 2*pc['O']
        del pc['O']
    if 'O' in bc.keys():
        bc['Fe2O3'] = bc['O']
        bc['FeO'] = bc['FeO'] - 2*bc['O']
        del bc['O']

#    pc_realMols = {}
    pc_mass = {}
    bc_mass = {}
    pc_tot = 0
    bc_tot = 0
    for e in pc.keys():
#        pc_realMols[e] = f*pc[e]
        pc_mass[e] = f*pc[e]*oxidemass[e]
        bc_mass[e] = bc[e]*oxidemass[e]
        pc_tot = pc_tot + pc_mass[e]
        bc_tot = bc_tot + bc_mass[e]

    X = pc_tot/bc_tot

    return X
