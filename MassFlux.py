#!/usr/bin/en3v python
"""
Author
------
Austin Smith (asmith155@alaska.edu)
"""

# Import supplemental modules.
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from scipy.stats import binned_statistic
import pickle as pkl


def compute_quad_areas_and_normals(x, y, z):
    # Ensure inputs are 2D numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    # Check dimensions
    if x.shape != y.shape or x.shape != z.shape:
        raise ValueError("x, y, and z arrays must have the same shape.")

    # Initialize lists to store face areas and normals
    face_areas = []
    face_normals = []

    # Loop over the grid to calculate areas and normals for each quadrilateral face
    for i in range(x.shape[0] - 1):
        for j in range(y.shape[1] - 1):
            # Define the corners of the quadrilateral
            p1 = np.array([x[i, j], y[i, j], z[i, j]])
            p2 = np.array([x[i + 1, j], y[i + 1, j], z[i + 1, j]])
            p3 = np.array([x[i + 1, j + 1], y[i + 1, j + 1], z[i + 1, j + 1]])
            p4 = np.array([x[i, j + 1], y[i, j + 1], z[i, j + 1]])

            # Debugging: print shapes
            #print(f"p1: {p1}, p2: {p2}, p3: {p3}, p4: {p4}")

            # Calculate the normal for one of the triangles (e.g., p1, p2, p4)
            edge1 = p2 - p1
            edge2 = p4 - p1
            
            #print(f"edge1: {edge1}, edge2: {edge2}")  # Debugging: print edges
            
            normal = np.cross(edge1, edge2)
            area = 0.5 * np.linalg.norm(normal)
            normal /= np.linalg.norm(normal)  # Normalize the normal vector

            # Calculate total area for the quadrilateral
            area_total = area * 2  # Since we have two triangles in the quadrilateral

            # Append the results
            face_areas.append(area_total)
            face_normals.append(normal)

    return face_areas, face_normals


with open("MFGridAndFluxes.pkl", 'rb') as file:
        SavedData = pkl.load(file)



               #number to kg      cm-3 to m^-3    km/s to m/s      Rj^2 to m^2
UnitConverter= 1.6726e-27 *       1e6 *           1e3*            (7.1492e7)**2

mfx= SavedData[3] #Density [cm-3] times Vx [km/s]
mfy= SavedData[4] #Density times Vy
mfz= SavedData[5] #Density times Vz

XX=SavedData[0] #[Rj]
YY=SavedData[1]
ZZ=SavedData[2]


MFperEgg=[]
AreaperEgg=[]
MaxR=[]
MinR=[]

for xcut in range(0,101):
    areas, normals = compute_quad_areas_and_normals(XX[xcut,:,:], YY[xcut,:,:],ZZ[xcut,:,:])
    MF=[]
    ilist=0
    for i in range(XX.shape[1] - 1):
        for j in range(YY.shape[2] - 1):
            MF.append((normals[ilist][0]*mfx[xcut,i,j]+normals[ilist][1]*mfy[xcut,i,j]+normals[ilist][2]*mfz[xcut,i,j])*areas[ilist]*UnitConverter)
            ilist = ilist +1
    
    MaxR.append(np.max(np.sqrt(XX[xcut,:,:]**2+ YY[xcut,:,:]**2+ ZZ[xcut,:,:]**2)))
    MinR.append(np.min(np.sqrt(XX[xcut,:,:]**2+ YY[xcut,:,:]**2+ ZZ[xcut,:,:]**2)))
    MFperEgg.append(sum(MF))
    AreaperEgg.append(sum(areas))

fig, ax1 = plt.subplots(figsize=(8, 6))
ax1.hlines(np.median(MFperEgg[10:40]),10,40,color='r',alpha=0.15)
ax1.plot(MFperEgg,c='k',label="Mass Flux",alpha=0.15)


ax2 = ax1.twinx() 
ax2.hlines(6,0,100,color="orange")
ax2.plot(MaxR,c='darkblue',label="Max Radius of Egg")
ax2.plot(MinR,c='lightblue',label="Min Radius of Egg")
ax1.set_ylim(-500,1500)
ax1.set_xlim(0,100)
ax2.set_ylim(0,50)
plt.legend()
ax1.grid()

ax2.grid(linestyle=":")
ax1.set_xlabel("Egg Cut")
ax1.set_ylabel("Mass Flux [kg/s]")
ax2.set_ylabel("Radius per Egg")
plt.savefig('MFperEgg.png', format='png')
#plt.show()
plt.close()










