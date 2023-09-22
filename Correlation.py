import numpy as np

FS_mass_bilan = np.loadtxt('Data\Fram_strait_mass_bilan.txt')
Index_A = np.loadtxt('Data\index_A.txt')


print(np.shape(FS_mass_bilan))
print(np.shape(Index_A))