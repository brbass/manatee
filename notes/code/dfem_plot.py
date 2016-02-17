import numpy as np
from matplotlib import pyplot as plt

phi = np.array([0.234451,
                0.234831,
                np.nan,
                0.237749,
                0.238122,
                np.nan,
                0.240995,
                0.241362,
                np.nan,
                0.244192,
                0.244554,
                np.nan,
                0.247342,
                0.247697,
                np.nan,
                0.250445,
                0.250795])

x = np.array([0.000000000000000000e+00,
              5.000000000000000104e-03,
              5.000000000000000104e-03,
              5.000000000000000104e-03,
              1.000000000000000021e-02,
              1.000000000000000021e-02,
              1.000000000000000021e-02,
              1.499999999999999944e-02,
              1.499999999999999944e-02,
              1.499999999999999944e-02,
              2.000000000000000042e-02,
              2.000000000000000042e-02,
              2.000000000000000042e-02,
              2.500000000000000139e-02,
              2.500000000000000139e-02,
              2.500000000000000139e-02,
              2.999999999999999889e-02])

z = np.array([0.000000000000000000e+00,
              5.000000000000000104e-03,
              1.000000000000000021e-02,
              1.499999999999999944e-02,
              2.000000000000000042e-02,
              2.500000000000000139e-02,
              2.999999999999999889e-02]) 

fig = plt.figure()
ax = fig.add_subplot(111, aspect=0.5)
ax.plot(x, phi)

plt.xlabel('r')
plt.ylabel('$\phi(r)$')
for val in z:
    plt.axvline(val, color='r')

plt.legend(('function', 'cell boundaries'), loc=4)

plt.savefig('dfem_plot.pdf', bbox_inches='tight')



