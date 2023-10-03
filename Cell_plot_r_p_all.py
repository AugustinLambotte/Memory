import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']

#montlhy
plt.figure(figsize = (15,7))
plt.suptitle('Correlation and p_value between sea ice melt\n and mean kinetic energy over the cells - monthly resolution',fontsize = 21)
for cell in Cell_name:    
    r = np.loadtxt('Data/'+cell+'/'+cell+'_r_monthly.txt')
    p = np.loadtxt('Data/'+cell+'/'+cell+'_p_monthly.txt')
    plt.subplot(121)
    plt.plot(np.arange(2011,2020),r,label = cell)
    plt.subplot(122)
    plt.plot(np.arange(2011,2020),p,label = cell)
plt.legend(fontsize = "15")
plt.grid()
plt.xlabel('year',fontdict = {'fontsize':20})
plt.ylabel('p',fontdict = {'fontsize':20})
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 14)
plt.subplot(121)
plt.legend(fontsize = "15")
plt.xlabel('year',fontdict = {'fontsize':20})
plt.ylabel('r',fontdict = {'fontsize':20})
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 14)
plt.grid()
plt.savefig('Plots\Correlation/4_cell/r_p_monthly_resolution.png')
plt.clf()
plt.close()
plt.figure(figsize = (15,7))

#daily
plt.suptitle('Correlation and p_value between sea ice melt\n and mean kinetic energy over the cells - daily resolution',fontsize = 21)
for cell in Cell_name:    
    r = np.loadtxt('Data/'+cell+'/'+cell+'_r_daily2.txt')
    p = np.loadtxt('Data/'+cell+'/'+cell+'_p_daily2.txt')
    if cell == 'Cell_B':
        plt.subplot(121)
        plt.plot(np.arange(2011,2020),r,label = cell,linewidth = 3.5)
        plt.subplot(122)
        plt.plot(np.arange(2011,2020),p,label = cell,linewidth = 3.5)
    else:
        plt.subplot(121)
        plt.plot(np.arange(2011,2020),r,label = cell)
        plt.subplot(122)
        plt.plot(np.arange(2011,2020),p,label = cell)
plt.legend(fontsize = "15")
plt.grid()
plt.xlabel('year',fontdict = {'fontsize':20})
plt.ylabel('p',fontdict = {'fontsize':20})
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 14)
plt.subplot(121)
plt.legend(fontsize = "15")
plt.xlabel('year',fontdict = {'fontsize':20})
plt.ylabel('r',fontdict = {'fontsize':20})
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 14)
plt.grid()
plt.savefig('Plots\Correlation/4_cell/r_p_daily_resolution.png')
