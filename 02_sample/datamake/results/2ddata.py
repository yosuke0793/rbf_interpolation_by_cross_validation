import numpy as np

data_full = np.loadtxt('dataset_full.txt', dtype=float)
print(data_full.shape)

data_2d = np.zeros((data_full.shape[0], 9))
data_2d[:, 0] = data_full[:, 0]
data_2d[:, 1] = data_full[:, 1]
data_2d[:, 2] = data_full[:, 3]
data_2d[:, 3] = data_full[:, 6]
data_2d[:, 4] = data_full[:, 7]
data_2d[:, 5] = data_full[:, 9]
data_2d[:, 6] = data_full[:, 12]
data_2d[:, 7] = data_full[:, 13]
data_2d[:, 8] = data_full[:, 15]

print(data_2d.shape)
np.savetxt('dataset.txt', data_2d, fmt='%16.6e')
