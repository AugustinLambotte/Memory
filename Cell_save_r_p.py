import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']

# Daily
for cell in Cell_name:
    Index = np.loadtxt('Data/'+cell+'/'+cell+'_index_daily.txt')
    MB = np.loadtxt('Data/'+cell+'/'+cell+'_MB_daily2.txt')
    r = [] #correlation coeficient 
    p = [] #p_value
    for year in range(len(Index)):
        if year == 1 or year == 5: #Bisextile year
            r.append(stats.pearsonr(Index[year,1:],MB[year,:])[0])
            p.append(stats.pearsonr(Index[year,1:],MB[year,:])[1])
        else:
            r.append(stats.pearsonr(Index[year,1:-1],MB[year,:-1])[0])
            p.append(stats.pearsonr(Index[year,1:-1],MB[year,:-1])[1])
    r = np.array(r)  
    p = np.array(p)  
    np.savetxt('Data/'+cell+'/'+cell+'_r_daily2.txt',r)
    np.savetxt('Data/'+cell+'/'+cell+'_p_daily2.txt',p)

# monthly
for cell in Cell_name:
    Index = np.loadtxt('Data/'+cell+'/'+cell+'_index_monthly.txt')
    MB = np.loadtxt('Data/'+cell+'/'+cell+'_MB_monthly.txt')
    r = [] #correlation coeficient 
    p = [] #p_value
    for year in range(len(Index)):
        r.append(stats.pearsonr(Index[year,:],MB[year,:])[0])
        p.append(stats.pearsonr(Index[year,:],MB[year,:])[1])
    r = np.array(r)  
    p = np.array(p)  
    np.savetxt('Data/'+cell+'/'+cell+'_r_monthly.txt',r)
    np.savetxt('Data/'+cell+'/'+cell+'_p_monthly.txt',p)

