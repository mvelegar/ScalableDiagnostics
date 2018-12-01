import numpy as np
import matplotlib.pyplot as plt
import timeit
import h5py


# Import color map
import cmocean

# import rsvd, rnmf, rspca from the ristretto pacakge
from ristretto.svd import rsvd
from ristretto.nmf import rnmf
from ristretto.pca import rspca



#******************************************************************************
# Read in Data
# Change path if needed....
#******************************************************************************    
f = h5py.File('data/O3.h5', 'r')
X = f['OH']

# This can take a while.... have a tea :)
t0 = timeit.default_timer()
X = np.asarray(X, order='F', dtype=np.float32)
X = X.T
print(timeit.default_timer()  - t0)


nLon=72
nLat=46 
nLev=30


# Sanity check, that we read in the correct data....
Xsnap = np.asarray(X[:,0])
Xsnap = Xsnap.reshape(nLon,nLat,nLev,  order='F') 
plt.imshow(Xsnap[:,::-1, 0].T, interpolation='bicubic', cmap=cmocean.cm.balance)
plt.colorbar()




#******************************************************************************
# i) Compute Randomized SVD
#******************************************************************************   
t0 = timeit.default_timer() 
U, S, Vt = rsvd(X, rank=250, oversample=20, n_subspace=3)
print(timeit.default_timer()  - t0)


# Plot Cum Engergy 
cum = np.cumsum(S**2) / np.linalg.norm(X)**2
plt.figure()
plt.plot(cum)
plt.show()


# Save decomposition...
hf = h5py.File('results/svd_OH.h5', 'w')
hf.create_dataset('U',   data=U, compression="gzip" )
hf.create_dataset('S',   data=S, compression="gzip")
hf.create_dataset('Vt',  data=Vt, compression="gzip")
hf.create_dataset('cum', data=cum, compression="gzip")
hf.close()



#******************************************************************************
# Plot the top 5 modes
#****************************************************************************** 
fig, axs = plt.subplots(2,5 , figsize=(20, 4), dpi=150, facecolor='w', edgecolor='k')
#fig.subplots_adjust(hspace = 1.5, wspace=1)
axs = axs.ravel()
for i in range(5):
    Xsnap = np.asarray(U[:,i])
    Xsnap = Xsnap.reshape(nLon, nLat, nLev,  order='F') 

    mode = Xsnap[:, ::-1, 0]


    im = axs[i].imshow(mode.T,  cmap=cmocean.cm.balance, interpolation='bicubic', vmin=-0.03, vmax=0.03)
    #axs[i].tick_params(axis='y', labelsize=14, color='white') 
    #axs[i].tick_params(axis='x', labelsize=14, color='white') 
    #axs[i].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')        
    #axs[i].grid('off')        
    fig.colorbar(im, ax=axs[i], orientation='vertical')


for i in range(5):
    time_dynamics = Vt[i,::360]
    axs[i+5].plot(time_dynamics, color='#252525')
    axs[i+5].tick_params(axis='y', labelsize=14) 
    axs[i+5].tick_params(axis='x', labelsize=14)         
    #axs[i+3].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')        
    axs[i+5].grid('off') 

fig.tight_layout()
plt.show()


del(U,S,Vt)




#******************************************************************************
# Compute Randomized NMF
#******************************************************************************   
t0 = timeit.default_timer()
W, H = rnmf(X, rank=8, oversample=200, n_subspace=3, maxiter=1000, verbose=True, tol=1e-5, init='nndsvd')
print("Total time: ", ( timeit.default_timer()  - t0 ))    


# Pull engergy 
C = np.sum(W**2, axis=0)
R = np.sum(H**2, axis=1)
S = C * R

# Plot Spectrum
plt.figure()
plt.plot(S)
plt.show()

# Reorder
idx = np.argsort(S)[::-1]
W = W[:,idx]
H = H[idx,:]
C = C[idx]
R = R[idx]
S = C**0.5 * R**0.5 

# Plot Cum Engergy 
cum = np.cumsum(S**2) / np.linalg.norm(X)**2
plt.figure()
plt.plot(cum)
plt.show()



hf = h5py.File('results/nmf_OH.h5', 'w')
hf.create_dataset('W',    data=W, compression="gzip" )
hf.create_dataset('H',    data=H, compression="gzip")
hf.create_dataset('S',    data=S, compression="gzip")
hf.create_dataset('cum',  data=cum, compression="gzip")

hf.close()




#******************************************************************************
# Plots
#****************************************************************************** 
fig, axs = plt.subplots(2,8 , figsize=(20, 4), dpi=150, facecolor='w', edgecolor='k')
#fig.subplots_adjust(hspace = 1.5, wspace=1)
axs = axs.ravel()

for i in range(8):
    Xsnap = np.asarray(W[:,i])
    Xsnap = Xsnap.reshape(nLon, nLat, nLev,  order='F') 
    
    #mode = Xsnap[:, ::-1, :].sum(axis=2)
    mode = Xsnap[:, ::-1, 0]
    #mode[mode < np.max(mode) * 0.5] = 0
    im = axs[i].imshow(mode.T,  cmap=cmocean.cm.amp, interpolation='bicubic')
    #axs[i].tick_params(axis='y', labelsize=14, color='white') 
    #axs[i].tick_params(axis='x', labelsize=14, color='white') 
    #axs[i].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')        
    #axs[i].grid('off')        
    fig.colorbar(im, ax=axs[i], orientation='vertical')


for i in range(8):
    time_dynamics = H[i,:]
    axs[i+8].plot(time_dynamics, color='#252525')
    axs[i+8].tick_params(axis='y', labelsize=14) 
    axs[i+8].tick_params(axis='x', labelsize=14)         
    #axs[i+3].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')        
    axs[i+8].grid('off') 

fig.tight_layout()
plt.show()







#******************************************************************************
# Compute Randomized Sparse PCA
#******************************************************************************   
from sklearn.preprocessing import scale


t0 = timeit.default_timer()
Bstar, Astar, eigvals, fit = rspca(scale(X.T, axis=0, with_mean=True, with_std=False, copy=True), n_components=25, 
                                   oversample=200, n_subspace=3, max_iter=200, 
                                   verbose=True, tol=1e-10, alpha = 1e-5, beta = 1e-12)
print("Total time: ", ( timeit.default_timer()  - t0 ))    

#plt.figure()
#plt.plot(fit)
#plt.show()


cum = np.cumsum(eigvals*(26208 - 1)) / np.linalg.norm(scale(X.T, axis=0, with_mean=True, with_std=False, copy=True))**2
plt.figure()
plt.plot(cum)
plt.show()




hf = h5py.File('results/spca_OH_less_sparse.h5', 'w')
hf.create_dataset('B',          data=Bstar, compression="gzip" )
hf.create_dataset('A',          data=Astar, compression="gzip")
hf.create_dataset('eigvals',    data=eigvals, compression="gzip")
hf.create_dataset('cum',        data=cum, compression="gzip")

hf.close()




#******************************************************************************
# Compute Randomized Sparse PCA
#******************************************************************************   
from sklearn.preprocessing import scale


t0 = timeit.default_timer()
Bstar, Astar, eigvals, fit = rspca(scale(X.T, axis=0, with_mean=True, with_std=False, copy=True), n_components=25, 
                                   oversample=200, n_subspace=3, max_iter=200, 
                                   verbose=True, tol=1e-10, alpha = 1e-4, beta = 1e-12)
print("Total time: ", ( timeit.default_timer()  - t0 ))    

#plt.figure()
#plt.plot(fit)
#plt.show()


cum = np.cumsum(eigvals*(26208 - 1)) / np.linalg.norm(scale(X.T, axis=0, with_mean=True, with_std=False, copy=True))**2
plt.figure()
plt.plot(cum)
plt.show()




hf = h5py.File('results/spca_OH.h5', 'w')
hf.create_dataset('B',          data=Bstar, compression="gzip" )
hf.create_dataset('A',          data=Astar, compression="gzip")
hf.create_dataset('eigvals',    data=eigvals, compression="gzip")
hf.create_dataset('cum',        data=cum, compression="gzip")

hf.close()







#******************************************************************************
# Plots
#****************************************************************************** 
fig, axs = plt.subplots(2,8 , figsize=(20, 4), dpi=150, facecolor='w', edgecolor='k')
#fig.subplots_adjust(hspace = 1.5, wspace=1)
axs = axs.ravel()
for i in range(8):
    Xsnap = np.asarray(Bstar[:,i])
    Xsnap = Xsnap.reshape(nLon, nLat, nLev,  order='F') 

    minmax = np.max(np.abs(Xsnap))
    im = axs[i].imshow(Xsnap[:, ::-1, 0].T,  cmap=cmocean.cm.balance, interpolation='bicubic', vmin=-minmax, vmax=minmax)
    #axs[i].tick_params(axis='y', labelsize=14, color='white') 
    #axs[i].tick_params(axis='x', labelsize=14, color='white') 
    #axs[i].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')        
    #axs[i].grid('off')        
    fig.colorbar(im, ax=axs[i], orientation='vertical')


Z = X[:,:].T.dot(Bstar)
Z = Z.T
for i in range(8):
    time_dynamics = Z[i,:]
    axs[i+8].plot(time_dynamics, color='#252525')
    axs[i+8].tick_params(axis='y', labelsize=14) 
    axs[i+8].tick_params(axis='x', labelsize=14)         
    #axs[i+3].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')        
    axs[i+8].grid('off') 

fig.tight_layout()
plt.show()










