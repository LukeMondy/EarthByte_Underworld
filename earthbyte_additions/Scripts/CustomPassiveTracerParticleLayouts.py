
# coding: utf-8

# # Using custom passive tracer particle layouts in Underworld
# 
# Underworld suffers from a poor algorithm for filling shapes with particles. This is a problem, since it can be useful to position passive tracers in shapes, such as spheres or walls, without a huge time penalty. 
# 
# To overcome this issue, this script allows users to define functions for particle layouts (the example given is for uniformly distributing points on a sphere). The points are then written to a h5 file, which Underworld can understand.
# 
# ## Point layout functions
# Any functions users define for points must return a numpy array of format:
# 
#     (x, y, (z), materialIndex)
# 
# The z axis is optional.
# 
# The resulting array of the above format can then be passed to the write_points_to_h5 function, which will handle the formatting automatically.
# 
# ## Using the particles in Underworld
# 
# In Underworld, you need to use a FileParticleLayout type. It only needs to take a 'filename' param. This param is the path and name of the data file you made when calling the write_points_to_h5 function - EXCEPT, it does not need the file extension. For example, "data.h5" would be simply "data".
# 
#     <struct name="particleLayout1"> 
#         <param name="Type"> FileParticleLayout </param>
#         <param name="filename"> <!-- the name of the file you made WITHOUT the .h5 extension --> </param>
#     </struct>
#     
# And Underworld will load the particles.
# 

# In[ ]:

try:
    import h5py
    import numpy as np
except ImportError as e:
    sys.exit("ERROR - You need to have h5py and numpy. Please install the one you are lacking. The computer says:\n%s" % e)


# In[ ]:

def write_points_to_h5(points, filename):
    size, n = points.shape
    
    pos = points[:,:n-1]
    mi = points[:,-1].astype(int)
    
    try:
        f = h5py.File(filename, "w")
        f.attrs['Swarm Particle Count'] = size

        dpos = f.create_dataset("Position", (size,n-1), dtype='f')
        dmi =  f.create_dataset("MaterialIndex", (size,), dtype='i')

        dpos[...] = pos
        dmi[...] = mi
    except Exception as e:
        sys.exit("\nERROR\nComputer says:\n%s\n" % e)
    finally:
        f.close()


# In[ ]:

def make_sphere_shell(num_points, centre_x, centre_y, centre_z, radius, materialIndex, dims):
    if 1 < dims <= 3:
    
        # Make points
        x = np.random.normal(size=(num_points, dims)) 
        x /= np.linalg.norm(x, axis=1)[:, np.newaxis]

        # Scale to radius
        x *= radius

        # Position the sphere
        x[:,0] += centre_x
        x[:,1] += centre_y
        if dims == 3:
            x[:,2] += centre_z
            
        # Tack on MaterialIndex
        mi = np.zeros((num_points,1))
        mi.fill(materialIndex)
        
        x = np.hstack((x, mi))

    else:
        sys.exit("\nERROR\nUnsupported number of dimensions. 2 or 3 is OK!\n")
            
    return x


def make_your_own_function_here(num_points, materialIndex, dims):
    points = np.zeros((num_points, dims + 1))
    
    # Put the details in here!
    
    points[:,-1] = materialIndex
    return points


# ## Example script 1
# The following shows an example script, using the functions defined above. 
# 
# It uses a loop to generate a number of spheres, and stores them all in the allPoints array. At the end, the allPoints array is written out using the write_points_to_h5 function

# In[ ]:

dims = 2                    # Make 2D spheres (circles)
num_spheres = 20            # How many spheres wanted
points_per_sphere = 2000    # How many points should each sphere have?
sphere_radii = 5000         # Radius

# np.linspace makes an evenly spaced set of numbers (40 in this case) between -185000 and 185000. 
xcentres = np.linspace(-185000, 185000, num_spheres)
# All the spheres will have the same Y value
ycentre = -140000

# Initialise the storage
allPoints = None

# For each x location from np.linspace
for count, xc in enumerate(xcentres):
    # make a sphere at that location
    points = make_sphere_shell(points_per_sphere, xc, ycentre, 0, sphere_radii, count, dims)
    
    # Stack all the loop outputs together
    if allPoints is None:
        allPoints = points
    else:
        allPoints = np.vstack((allPoints, points))

filename = "2Ddata.h5"
write_points_to_h5(allPoints, filename)


# In[ ]:

from matplotlib import pyplot as plt
plt.figure(figsize=(20,20))
plt.scatter(allPoints[:,0], allPoints[:,1], c=allPoints[:,2])

xmin, xmax = min(allPoints[:,0]), max(allPoints[:,0])
ymin, ymax = min(allPoints[:,1]), max(allPoints[:,1])

xbuff = (xmax - xmin) / 10.0
ybuff = (ymax - ymin) / 5.0

plt.xlim(xmin - xbuff, xmax + xbuff)
plt.ylim(ymin - ybuff, ymax + ybuff)

plt.axes().set_aspect('equal')
plt.savefig("2D_plot.png")


# ## Example script 2
# The following shows another example script, now using 3D spheres, in a 3D layout.

# In[ ]:

dims = 3                    # Make 3D spheres (circles)
x_num_spheres = 10          # How many spheres wanted in X direction
y_num_spheres =  4          # How many spheres wanted in Y direction (Y is up and down)
z_num_spheres = 10          # How many spheres wanted in Z direction

points_per_sphere = 1000    # How many points should each sphere have?
sphere_radii = 5000         # Radius

# np.linspace makes an evenly spaced set of numbers (10 in this case) between -185000 and 185000. 
xcentres = np.linspace(-185000, 185000, x_num_spheres)
ycentres = np.linspace(-40000, 0, y_num_spheres)
zcentres = np.linspace(-185000, 185000, z_num_spheres)

# Initialise the storage
allPoints = None
count = 0

# For each set of centres, we need to make a circle - that is x_num_spheres * y_num_spheres * z_num_spheres
# which equals 400 spheres in this case. 
# Any order of the for loops should work, but typically the outer loop is for vertical axes.
for yc in ycentres:
    for zc in zcentres:
        for xc in xcentres:
            # make a sphere at that location
            points = make_sphere_shell(points_per_sphere, xc, yc, zc, sphere_radii, count, dims)
            
            # Stack all the loop outputs together
            if allPoints is None:
                allPoints = points
            else:
                allPoints = np.vstack((allPoints, points))
                
            count += 1

filename = "3Ddata.h5"
write_points_to_h5(allPoints, filename)


# ### Write a XMF file to look at the data in Paraview
# Paraview can visualise the 3D point data we made once it knows the general structure of the data.
# 
# The following XML code gives Paraview enough to work with. Open paraview, and load the 3D_paraview.xmf file to see the particles you made.

# In[ ]:

xmf = """<?xml version="1.0" ?>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
    <Domain>
        <Grid Name="materialSwarm" GridType="Collection">
            <Time Value="0" />
            <Grid Name="spheres">
                <Topology Type="POLYVERTEX" NodesPerElement="{numPoints}"> </Topology>
                <Geometry Type="XYZ">
                    <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="{numPoints} 3">{filename}:/Position</DataItem>
                </Geometry>

                <Attribute Type="Scalar" Center="Node" Name="MaterialIndex">
                    <DataItem Format="HDF" NumberType="Int" Dimensions="{numPoints} 1">{filename}:/MaterialIndex</DataItem>
                </Attribute>
            </Grid>
        </Grid>
    </Domain>
</Xdmf>
""".format(numPoints=allPoints.shape[0], filename=filename)

try:
    with open("3D_paraview.xmf", 'w') as f:
        f.write(xmf)
except Exception as e:
    sys.exit("ERROR with writing 3D_paraview.xmf. Computer says:\n%s" % e)

