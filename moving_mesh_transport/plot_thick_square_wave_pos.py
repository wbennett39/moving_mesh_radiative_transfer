import matplotlib.pyplot as plt
import numpy as np


tpnts = [1, 5, 10, 25]
wave_loc = [400.4, 400.1, 400, 400] 

plt.plot(tpnts, wave_loc, '-o')
# plt.plot(tpnts, 2/np.sqrt(tpnts), 'k-')
plt.show()