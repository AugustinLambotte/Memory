import numpy as np
import matplotlib.pyplot as plt

FS_mass_bilan = np.loadtxt('Data\Fram_strait_mass_bilan.txt')
Index_A = np.loadtxt('Data\index_A.txt')
Ricker = np.loadtxt('Data\Ricker.txt')
M1 = np.loadtxt('Data\M1.txt')
M2 = np.loadtxt('Data\M2.txt')
r_coef = np.zeros(len(FS_mass_bilan))
r_coef_M2 = np.zeros(np.shape(M2)[0])
month_names = ["jan", "feb", "mar", "april", "may", "june", "july", "aug", "sept", "oct","nov", "dec"]

for i in range(len(FS_mass_bilan)):
    r_coef[i] = np.corrcoef(FS_mass_bilan[i],Index_A[i])[0,1]

for i in range(np.shape(M2)[0]):
    r_coef_M2[i] = np.corrcoef(M2[i],Index_A[i])[0,1]

fig = plt.figure(figsize = (12,7))
fig.suptitle('Correlation coefficient betwwen Index A and Fram Strait mass volume data',fontsize = 20)
ax = plt.axes()
ax.plot([year for year in range(2011,2020)],r_coef,label = 'myself')
ax.plot([year for year in range(2011,2017)],r_coef_M2,label = 'M2')
ax.grid()
ax.legend(fontsize = 15)
ax.tick_params(axis='y',labelsize = 20)
ax.tick_params(axis='x',labelsize = 20)
ax.set_xlabel('year',fontsize = 20)
ax.set_ylabel('r',fontsize = 20)

plt.savefig('Plots/Correlation/r_evolution.png')
plt.clf()

quit()

for year in range(2011,2020):
    #Param
    fig.suptitle(f'{year}',fontsize = 25)
    ax1 = plt.axes()
    ax2 = ax1.twinx()
    ax1.grid()

    bbox = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)
    ax1.text(0.95, 0.07, f"r = {round(r_coef[year-2011],2)}", fontsize=15, bbox=bbox,
            transform=ax1.transAxes, horizontalalignment='right')

    ax1.set_ylabel('km^3',color = 'red',fontsize = 20)
    ax2.set_ylabel('J/kg',color = 'blue',fontsize = 20)
    ax1.tick_params(axis='y',labelcolor = 'red',labelsize = 20)
    ax2.tick_params(axis='y',labelcolor = 'blue',labelsize = 20)
    ax1.tick_params(axis='x',labelsize = 20)
    ax1.set_xlabel('month',fontdict={'fontsize':20})
    #Plot
    lns1 = ax1.plot([month for month in range(1,13)],FS_mass_bilan[year-2011,:], label = "FS_mass_bilan",color = 'red')
    lns2 = ax2.plot([month for month in range(1,13)],Index_A[year-2011,:],label = 'Index_A', color = 'blue')

    #legend
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, fontsize = 15)

    plt.savefig(f'Plots/Correlation/FS-A/{year}_correlation.png')
    plt.clf()

