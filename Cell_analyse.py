import numpy as np
import matplotlib.pyplot as plt

Cell_name = ['Cell_A','Cell_B','Cell_C','Cell_D']

# Daily
for cell in Cell_name:
    Index = np.loadtxt('Data/'+cell+'/'+cell+'_index_daily.txt')
    MB = np.loadtxt('Data/'+cell+'/'+cell+'_MB_daily.txt')
    r = [] #correlation coeficient 
    for year in range(len(Index)):
        if year == 1 or year == 5: #Bisextile year
            r.append(np.corrcoef(Index[year,1:],MB[year,:])[0,1])
        else:
            r.append(np.corrcoef(Index[year,1:-1],MB[year,:-1])[0,1])
    plt.plot(np.arange(2011,2020),r,label = cell)
plt.legend()
plt.title('daily')

plt.grid()
plt.show()

#Monthly
for cell in Cell_name:
    Index = np.loadtxt('Data/'+cell+'/'+cell+'_index_monthly.txt')
    MB = np.loadtxt('Data/'+cell+'/'+cell+'_MB_monthly.txt')
    r = [] #correlation coeficient 
    for year in range(len(Index)):
        r.append(np.corrcoef(Index[year,:],MB[year,:])[0,1])
    plt.plot(np.arange(2011,2020),r,label = cell)
plt.legend()
plt.title('Monthly')
plt.grid()
plt.show()