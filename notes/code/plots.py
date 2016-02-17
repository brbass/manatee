import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

text = ""

u235 = np.loadtxt("u235.out", skiprows=2)
ud2o = np.loadtxt("ud20.out", skiprows=2)
uran = np.loadtxt("uranium.out",skiprows=2)

#########
# U-235 #
#########

x = np.unique(u235[:,0])
o = np.unique(u235[:,1])

num_x = len(x)
num_o = len(o)

it = np.zeros([num_x, num_o])
t = np.zeros([num_x, num_o])
k = np.zeros([num_x, num_o])

index = 0
index_x = 0
for xval in x:
    index_o = 0
    for oval in o:
        it[index_x, index_o] = u235[index, 2]
        t[index_x, index_o] = u235[index, 3]
        k[index_x, index_o] = abs(u235[index, 4] - 1)
        
        index += 1
        index_o += 1
    index_x += 1

plt.figure()
plt.pcolormesh(o, x, k, norm=mpl.colors.LogNorm(vmin=k.min().min(), vmax=k.max().max()))
#plt.title('U-235, Abs error in k')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("k_u235{0}.pdf".format(text),bbox_inches='tight')

plt.figure()
plt.pcolormesh(o, x, t, norm=mpl.colors.LogNorm(vmin=t.min().min(), vmax=t.max().max()))
#plt.title('U-235, Solve time (sec)')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("t_u235{0}.pdf".format(text),bbox_inches='tight')

plt.figure()
plt.pcolormesh(o, x, it)#, norm=mpl.colors.LogNorm(vmin=it.min().min(), vmax=it.max().max()))
#plt.title('U-235, Number of iterations')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("it_u235{0}.pdf".format(text),bbox_inches='tight')

#########
# U-D2O #
#########

x = np.unique(ud2o[:,0])
o = np.unique(ud2o[:,1])

num_x = len(x)
num_o = len(o)

it = np.zeros([num_x, num_o])
t = np.zeros([num_x, num_o])
k = np.zeros([num_x, num_o])

index = 0
index_x = 0
for xval in x:
    index_o = 0
    for oval in o:
        it[index_x, index_o] = ud2o[index, 2]
        t[index_x, index_o] = ud2o[index, 3]
        k[index_x, index_o] = abs(ud2o[index, 4] - 1)
        
        index += 1
        index_o += 1
    index_x += 1

plt.figure()
plt.pcolormesh(o, x, k, norm=mpl.colors.LogNorm(vmin=k.min().min(), vmax=k.max().max()))
#plt.title('U-D2O, Abs error in k')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("k_ud2o{0}.pdf".format(text),bbox_inches='tight')

plt.figure()
plt.pcolormesh(o, x, t, norm=mpl.colors.LogNorm(vmin=t.min().min(), vmax=t.max().max()))
#plt.title('U-D2O, Solve time (sec)')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("t_ud2o{0}.pdf".format(text),bbox_inches='tight')

plt.figure()
plt.pcolormesh(o, x, it)#, norm=mpl.colors.LogNorm(vmin=it.min().min(), vmax=it.max().max()))
#plt.title('U-D2O, Number of iterations')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("it_ud2o{0}.pdf".format(text),bbox_inches='tight')

###########
# Uranium #
###########

x = np.unique(uran[:,0])
o = np.unique(uran[:,1])

num_x = len(x)
num_o = len(o)

it = np.zeros([num_x, num_o])
t = np.zeros([num_x, num_o])
k = np.zeros([num_x, num_o])

index = 0
index_x = 0
for xval in x:
    index_o = 0
    for oval in o:
        it[index_x, index_o] = uran[index, 2]
        t[index_x, index_o] = uran[index, 3]
        k[index_x, index_o] = abs(uran[index, 4] - 1)
        
        index += 1
        index_o += 1
    index_x += 1

plt.figure()
plt.pcolormesh(o, x, k, norm=mpl.colors.LogNorm(vmin=k.min().min(), vmax=k.max().max()))
#plt.title('Uranium, Abs error in k')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("k_uran{0}.pdf".format(text),bbox_inches='tight')

plt.figure()
plt.pcolormesh(o, x, t, norm=mpl.colors.LogNorm(vmin=t.min().min(), vmax=t.max().max()))
#plt.title('Uranium, Solve time (sec)')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("t_uran{0}.pdf".format(text),bbox_inches='tight')

plt.figure()
plt.pcolormesh(o, x, it)#, norm=mpl.colors.LogNorm(vmin=it.min().min(), vmax=it.max().max()))
#plt.title('Uranium, Number of iterations')
plt.xlabel('Number of ordinates')
plt.ylabel('Number of cells')
cb = plt.colorbar()
plt.savefig("it_uran{0}.pdf".format(text),bbox_inches='tight')

