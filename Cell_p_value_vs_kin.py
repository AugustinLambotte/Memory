import numpy as np
import matplotlib.pyplot as plt

Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']


#Total kin_ene
for cell in Cell_name:
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('r',color = 'red')
    ax2.set_ylabel('kin_ene',color = 'blue')
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')
    r = np.loadtxt(f'Data\{cell}\{cell}_r_daily2.txt')
    kin_ene = np.loadtxt(f'Data\{cell}\{cell}_index_daily.txt')
    kin_ene = np.nansum(kin_ene,axis=1)
    ax1.plot(np.arange(2011,2020),r,color = 'red')
    #ax2.plot(np.arange(2011,2020),melting_brut,color = 'blue',linestyle = 'dashed')
    ax2.plot(np.arange(2011,2020),kin_ene,color = 'blue')
    
    fig.suptitle(cell)
    ax1.grid()
    plt.show()
    plt.close()

#mean kin_ene
for cell in Cell_name:
    fig = plt.figure()
    ax1 = plt.axes()
    ax2 = plt.twinx()
    ax1.set_xlabel('year')
    ax1.set_ylabel('r',color = 'red')
    ax2.set_ylabel('kin_ene',color = 'blue')
    ax1.tick_params(axis = 'y', color = 'red', labelcolor = 'red')
    ax2.tick_params(axis = 'y', color = 'blue', labelcolor = 'blue')
    r = np.loadtxt(f'Data\{cell}\{cell}_r_daily2.txt')
    kin_ene = np.loadtxt(f'Data\{cell}\{cell}_index_daily.txt')
    kin_ene = np.nanmean(kin_ene,axis=1)
    ax1.plot(np.arange(2011,2020),r,color = 'red')
    #ax2.plot(np.arange(2011,2020),melting_brut,color = 'blue',linestyle = 'dashed')
    ax2.plot(np.arange(2011,2020),kin_ene,color = 'blue')
    
    fig.suptitle(cell)
    ax1.grid()
    plt.show()
    plt.close()
