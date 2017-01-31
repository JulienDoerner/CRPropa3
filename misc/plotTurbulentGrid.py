# CRPropa example
# Evaluates and plots the turbulent fields generated by CRPropa

from crpropa import *
from pylab import *


# Create turbulent field
origin = Vector3d(0,0,0)
n = 128
spacing = 1
lMin, lMax = 2, 32
Brms = 1
alpha = -11./3.
seed = 42
Lc = turbulentCorrelationLength(lMin, lMax, alpha)
vGrid = VectorGrid(origin, n, spacing)
initTurbulence(vGrid, Brms, lMin, lMax, alpha, seed)

# Copy field grid to array(s)
Bx, By, Bz = zeros((3, n, n, n))
for ix in range(n):
    for iy in range(n):
        for iz in range(n):
            b = vGrid.get(ix, iy, iz)
            Bx[ix, iy, iz] = b.x
            By[ix, iy, iz] = b.y
            Bz[ix, iy, iz] = b.z

# Plot slice in position space
i = 0  #n/2
sxy = (Bx[:,:,i]**2 + By[:,:,i]**2 + Bz[:,:,i]**2)**.5
sxz = (Bx[:,i,:]**2 + By[:,i,:]**2 + Bz[:,i,:]**2)**.5
syz = (Bx[i,:,:]**2 + By[i,:,:]**2 + Bz[i,:,:]**2)**.5

figure(figsize=(8,6))
im = imshow(sxy, origin='lower', extent=[0,n,0,n], vmin=0, vmax=3)
cbar = colorbar(im)
cbar.set_label('$B/B_{rms}$')
xlabel('$x$')
ylabel('$y$')
text(0.8, 1.05, '$z=%i$'%(n/2), transform=gca().transAxes)
savefig('TurbulentGrid_slicePositionSpace_xy.png', bbox_inches='tight')

figure(figsize=(8,6))
im = imshow(sxz, origin='lower', extent=[0,n,0,n], vmin=0, vmax=3)
cbar = colorbar(im)
cbar.set_label('$B/B_{rms}$')
xlabel('$x$')
ylabel('$z$')
text(0.8, 1.05, '$z=%i$'%(n/2), transform=gca().transAxes)
savefig('TurbulentGrid_slicePositionSpace_xz.png', bbox_inches='tight')

figure(figsize=(8,6))
im = imshow(syz, origin='lower', extent=[0,n,0,n], vmin=0, vmax=3)
cbar = colorbar(im)
cbar.set_label('$B/B_{rms}$')
xlabel('$y$')
ylabel('$z$')
text(0.8, 1.05, '$z=%i$'%(n/2), transform=gca().transAxes)
savefig('TurbulentGrid_slicePositionSpace_yz.png', bbox_inches='tight')

"""
# Plot slice in configuration space
figure()
Bkx = fftshift(fftn(Bx))
Bky = fftshift(fftn(By))
Bkz = fftshift(fftn(Bz))
Bk = (((Bkx*Bkx.conjugate() + Bky*Bky.conjugate() + Bkz*Bkz.conjugate()).real)**.5)
im = imshow(log10(Bk[:,:,n/2]), origin='lower', vmin=0)
cbar = colorbar(im)
cbar.set_label(r'$\log_{10}(B/B_{rms})$')
xlabel('$k_x$')
ylabel('$k_y$')
k = fftshift(fftfreq(n))
idx = arange(0, n, n/4)
xticks(idx, k[idx])
yticks(idx, k[idx])
xlim(0,n)
ylim(0,n)
text(0.8, 1.05, '$k_z=%.2f$'%k[n/2], transform=gca().transAxes)
savefig('TurbulentGrid_sliceConfigurationSpace.png', bbox_inches='tight')

# Histogram of field strengths
figure()
Bx.resize(n**3)
By.resize(n**3)
Bz.resize(n**3)
hist(log10(Bx), bins=40, range=(-3,3), histtype='step', normed=True, label='$B_x$', linewidth=2)
hist(log10(By), bins=40, range=(-3,3), histtype='step', normed=True, label='$B_y$', linewidth=2)
hist(log10(Bz), bins=40, range=(-3,3), histtype='step', normed=True, label='$B_z$', linewidth=2)
legend()
grid()
xlabel('$\log_{10}(B/B_{RMS}$)')
ylabel('Rel. Frequency')
brms = (mean( Bx**2 + By**2 + Bz**2 ))**.5
bmean = abs(mean(Bx + By + Bz))
text(0.95, 0.7, 'RMS = %.2f\nMean = %.2f'%(brms, bmean), ha='right', va='top', transform=gca().transAxes)
savefig('TurbulentGrid_amplitude.png', bbox_inches='tight')

show()
"""
