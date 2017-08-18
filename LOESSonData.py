# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 00:37:34 2017

@author: The Gags
"""
'''
Hi. 
Find the User Input Area, to change the Smoothness Factor or change the
input worksheet.
You will find the Smooth Normal points in x_nor and res_normal. x_nor contains the 
x coordinates. res_normal contains the y coordinates.
You will find the Smooth Infected points in x_inf and res_inf. x_nor contains the 
x coordinates. res_inf contains the y coordinates.

Also notice, in the "Computation" part, you can change the for loop's range to
include more genes. Right now it is set to range(0,1) to run ONLy for the first
gene.
'''




#START LOESS ******************************************************************

#START imports-----------------------------------------------------------------
import numpy as np
import statsmodels.api as sm
import pylab as plt
import xlrd as xl
import itertools as it
import math
#END imports-------------------------------------------------------------------

#START LOESS and Swap Functions------------------------------------------------
def loess ( x, w, l, xw, yw):
    xnew = []
    ynew = []
    left = 0
    right = 0
    for i in range(l):
        left = int(i - w / 2)
        #print("window:",window)
        
        if(left < 0):
            left = 0
        
        right = w + left
        #print("right:",right)
        if(right > l):
            left = l - w
            right = l
        
        #print("left:",left," right:",right)
            
        xnew = np.array(xw[left:right])          
        ynew = np.array(yw[left:right])  
    xdist = abs(x - xnew)
    #print("xdist:",xdist)
    xmax = max(xdist)
    xscaled = xdist / xmax
    wt = pow((1-pow(xscaled,3)),3)
    xnew = sm.add_constant(xnew)
    model = sm.WLS(ynew, xnew, wt)
    tup = model.fit_regularized().params      
    return(tup[1]*x + tup[0])

def swap( A, x, y ):
  tmp = A[x]
  A[x] = A[y]
  A[y] = tmp
#END LOESS and Swap Functions--------------------------------------------------

#START User Input Area---------------------------------------------------------  
alpha = 0.3  
wb = xl.open_workbook('DSSample.xlsx')
#col 0-14 is normal
#col 15-53 is infected 
#END User Input Area-----------------------------------------------------------

#START reading DataSet---------------------------------------------------------
for s in wb.sheets():
    #print 'Sheet:',s.name
    values = []
    sample_id = []
    gene_id = []
    
    for row in range(1,s.nrows):
        col_value = []
        for col in range(2,s.ncols):
            value  = (s.cell(row,col).value)
            try : value = float(value)
            except : pass
            col_value.append(value)
        values.append(col_value)
    
    for col in range(2,s.ncols):
        sample_id.append(str(s.cell(0,col).value))
    for row in range(1,s.nrows):
        gene_id.append(str(s.cell(row,1).value))
        
#END Reading Data Set----------------------------------------------------------
        
        
#START Computation ------------------------------------------------------------        
for gene in range(0,1):    
    
    #Normal Data
    x_nor = []
    y_nor = []
    nor_length = 0 
    
    for pair in it.combinations(values[gene][0:14], 2):
        x_nor.append((pair[0] + pair[1])/2)
        y_nor.append(math.sqrt(pair[0] * pair[1]))
        nor_length = nor_length + 1
    
    for i in range(nor_length):                 #Sorting
        for j in range(nor_length - 1 - i):
            if(x_nor[j]>x_nor[j+1]):
                swap(x_nor, j, j+1)
                swap(y_nor, j, j+1)
    
    print("Normal length:", nor_length)   
    window = int((nor_length*alpha))
    res_normal = []
    
    for i in range(nor_length):
        res_normal.append(loess(x_nor[i], window, nor_length, x_nor, y_nor))
        
    
    #Infected Data
    x_inf = []
    y_inf = []
    inf_length = 0
    
    for pair in it.combinations(values[0][15:53], 2):
        x_inf.append((pair[0] + pair[1])/2)
        y_inf.append(math.sqrt(pair[0] * pair[1]))
        inf_length = inf_length + 1
        
    for i in range(inf_length):                 #Sorting
        for j in range(inf_length - 1 - i):
            if(x_inf[j]>x_inf[j+1]):
                swap(x_inf, j, j+1)
                swap(y_inf, j, j+1)
    
    print("Infected length:", inf_length)
    window = int((inf_length*alpha))
    res_inf = []
    
    for i in range(inf_length):
        res_inf.append(loess(x_inf[i], window, inf_length, x_inf, y_inf))
#END Computation---------------------------------------------------------------    
 
#START Graph Plotting----------------------------------------------------------   
    
    plt.title(gene_id[gene])
    plt.plot(x_nor,y_nor,"+", label = "Normal Data")
    #plt.title("Normal")
    plt.plot(x_nor,res_normal, label = "Normal Trend")
    #plt.xscale('log')
    #plt.show()
    
    plt.plot(x_inf,y_inf,'.', label = "Infected Data")
    #plt.title("Infected")
    plt.plot(x_inf,res_inf, label = "Infected Trend")
    #plt.xscale('log')
    plt.legend()
    plt.show()
    
    #Commented out because you only need smoothed out points and noy the graphs
    
#END Graph Plotting------------------------------------------------------------

#END LOESS ********************************************************************




    



