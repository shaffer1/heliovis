import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pymeshlab
from scipy.spatial import Delaunay

# Download data set from plotly repo
#pts = np.loadtxt(np.DataSource().open('https://raw.githubusercontent.com/shaffer1/heliovis/main/minb.txt'))
#x, y, z = pts.T

#fig = go.Figure(data=[go.Mesh3d(x=x, y=y, z=z, color='lightpink', opacity=0.50)])
#fig.show()
import plotly.graph_objects as go

class timeStep:
  def __init__(self, timeID, minBs, start, stop,name):
    self.timeID = timeID
    self.name = name + str(timeID)
    self.ms = pymeshlab.MeshSet()
    self.minBs = minBs[start:stop]
    self.x = self.minBs[0:,0]
    self.y = self.minBs[0:,1]
    self.z = self.minBs[0:,2]

    self.minValx = np.min(self.x)
    self.maxValx = np.max(self.x)
    self.meanValx = np.mean(self.x)
    self.varValx = np.var(self.x)

    self.minValy = np.min(self.y)
    self.maxValy = np.max(self.y)
    self.meanValy = np.mean(self.y)
    self.varValy = np.var(self.y)

    self.minValz = np.min(self.z)
    self.maxValz = np.max(self.z)
    self.meanValz = np.mean(self.z)
    self.varValz = np.var(self.z)

    self.ssdz= self.sumSqFromZ0()



  def print(self):
      print("Timestep ", self.timeID," Num points: ", self.minBs.shape[0])
      print(self.minValz, " ", self.maxValz, " ", self.meanValz, " ", self.varValz)
      print("Sum of sqaured distances from Z0: ", self.ssdz)

  def triangulate2D(self):
      pts2D= self.minBs[0:,0:2]
      #print(pts2D)
      self.faces = Delaunay(pts2D).simplices
      self.mesh = pymeshlab.Mesh(self.minBs, self.faces)
      self.ms.add_mesh(self.mesh,self.name)
      return None
  
  def sumSqFromZ0(self):
      ssdz = self.minBs[0:,2]*self.minBs[0:,2]
      return np.sum(ssdz)
  
  def sumSqFromLSBF(self):
      return None
  
  def saveToOBJ(self,filename):
      if (self.mesh == None):
          self.triangulate2d()
      self.ms.save_current_mesh(filename)
      
      
  
  
   

def find_lsbf_plane(coords):
    # barycenter of the points
    # compute centered coordinates
    G = coords.sum(axis=0) / coords.shape[0]

    # run SVD
    u, s, vh = np.linalg.svd(coords - G)

    # unitary normal vector
    u_norm = vh[2, :]
    return u_norm

def find_dist_z0(coords):
    total_z = 0;
    i=0
    dists = np.zeros(coords.shape[0])
    for pt in coords:
        dists[i]=np.abs(pt[2])
        i=i+1
        total_z+=np.abs(pt[2])
    print("Total number points", coords.shape[0])
    print("Average distance: ", np.average(dists))
    print("Standard deviation: ", np.std(dists))
    print("Validation ", total_z/coords.shape[0])

def create_timesteps(coords, numsteps, stepSet):
    for i in range(numsteps):
        print("Making step", i)
        stepSet.append(timeStep(i,coords, i*480, (i+1)*480, "HeidiFieldLine"))

    
stepSet = []
print("Start it up")
# Text file data converted to integer data type
File_data = np.loadtxt("HeidiFieldLineXYZ2.txt", dtype=float,comments="#", skiprows=0, usecols=(2,3,4))
#print(File_data.shape)

field_lines = np.array_split(File_data,File_data.shape[0]/5)
#print(field_lines)
print("Processing ", len(field_lines), " field lines")
minb_pts=np.zeros((len(field_lines),3))
line_num=0
for line in field_lines:
    minb_pts[line_num] = line[2]
    line_num=line_num+1
#print(minb_pts)
find_dist_z0(minb_pts)

create_timesteps(minb_pts,540,stepSet)
tsIndex=-1
tsCurIndex=0
maxVar=0
maxssdz=0
tsSSDz=stepSet[0]
tsMaxVar = stepSet[0] 
for t in stepSet:
    if t.varValz > maxVar:
        tsIndex=tsCurIndex
        maxVar=t.varValz
        tsMaxVar=t
    if t.ssdz> maxssdz:
        maxssdz=t.ssdz
        tsSSDz=t
    tsCurIndex=tsCurIndex+1

print("step ", tsIndex, " has variance ", maxVar)
tsMaxVar.print()
tsMaxVar.sumSqFromZ0()
tsSSDz.print()
tsSSDz.triangulate2D()
tsSSDz.saveToOBJ("tstest.obj")



#print(field_lines)

"""
fig = go.Figure(data=go.Scatter3d())

i=0
for line in field_lines:
    #print("Line ", i)
    if (i==0):
        split_line = np.transpose(line)
    #print(line)
        xl = split_line[0] 
        yl = split_line[1]
        zl = split_line[2]
        fig= go.Figure(data=go.Scatter3d(x=xl, y=yl,z=zl, mode='lines'))
    else:
        split_line = np.transpose(line)
        xl = split_line[0] 
        yl = split_line[1]
        zl = split_line[2]
        fig.add_trace(go.Scatter3d(x=xl, y=yl,z=zl, mode='lines'))
    i=i+1
    #if (i==2):
     #       break

 fig = go.Figure(data=go.Cone(
    x=[1, 2, 3],
    y=[1, 2, 3],
    z=[1, 2, 3],
    u=[1, 0, 0],
    v=[0, 3, 0],
    w=[0, 0, 2],
    sizemode="absolute",
    sizeref=2,
    anchor="tip")) 

fig.update_layout(
      scene=dict(domain_x=[0, 1],
                 camera_eye=dict(x=-1.57, y=1.36, z=0.58)))

fig.show()
"""