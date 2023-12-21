import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pymeshlab
from scipy.spatial import Delaunay
import vtkmodules.all as vtk
import vtkmodules.util.numpy_support as numpy_support
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.numpy_interface import dataset_adapter as dsa
import csv

# Download data set from plotly repo
#pts = np.loadtxt(np.DataSource().open('https://raw.githubusercontent.com/shaffer1/heliovis/main/minb.txt'))
#x, y, z = pts.T

#fig = go.Figure(data=[go.Mesh3d(x=x, y=y, z=z, color='lightpink', opacity=0.50)])
#fig.show()
import plotly.graph_objects as go

#importing the os module
import os

#to get the current working directory
directory = os.getcwd()

print(directory)

class timeStep:
  def __init__(self, timeID, minBs, start, stop,name):
    self.timeID = timeID
    self.name = name + str(timeID)
    self.ms = pymeshlab.MeshSet()
    self.minBs = minBs[start:stop]
    self.x = self.minBs[0:,0]
    self.y = self.minBs[0:,1]
    self.z = self.minBs[0:,2]
    self.numVertices = self.z.shape[0]

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
    self.triangulate2D()
    self.numFaces = self.faces.shape[0]
    #print("N V ", self.numVertices, " N F ", self.numFaces)

  def get_max_z(self):
      return np.max(np.absolute(self.z))

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
  
  def saveToFile(self,filename):
      self.ms.save_current_mesh(filename)
      print("Writing out ", filename)

  def exportSurfaceDataToVTK(self,filename):
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename+".vtk")
    writer.SetInputData(self.my_vtk_polydata)
    writer.Write()
  
  def exportToVTK(self, filename):

    my_vtk_dataset = vtk.vtkUnstructuredGrid()

    points = vtk.vtkPoints()
    for id in range(self.numVertices):
        points.InsertPoint(id, [self.x[id], self.y[id], self.z[id]])
    my_vtk_dataset.SetPoints(points)


    my_vtk_dataset.Allocate(self.numFaces)
    for id in range(self.numFaces):
        point_ids = [self.faces[id,0], self.faces[id,1],self.faces[id,2]]
        my_vtk_dataset.InsertNextCell(vtk.VTK_TRIANGLE, 3, point_ids)

    m = self.ms.current_mesh()
    density = m.vertex_custom_scalar_attribute_array("Density")
    pressure = m.vertex_custom_scalar_attribute_array("Pressure") 

    data_type = vtk.VTK_DOUBLE
    vtk_density = numpy_support.numpy_to_vtk(num_array=density, deep=True, array_type=data_type)
    vtk_density.SetName("Density")
    my_vtk_dataset.GetPointData().AddArray(vtk_density)

    vtk_pressure = numpy_support.numpy_to_vtk(num_array=pressure, deep=True, array_type=data_type)
    vtk_pressure.SetName("Pressure")
    my_vtk_dataset.GetPointData().AddArray(vtk_pressure)
    
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename+".vtu")
    writer.SetInputData(my_vtk_dataset)
    writer.Write()
    print("Wrote ", filename+".vtu")
    
    return None
  
  def computeCurvature(self):
    print("Computing Curvature")
    self.my_vtk_polydata = vtk.vtkPolyData()
    triangles = vtk.vtkCellArray()
    points = vtk.vtkPoints()

    for id in range(self.numVertices):
        points.InsertPoint(id, [self.x[id], self.y[id], self.z[id]])
    self.my_vtk_polydata.SetPoints(points)

    for id in range(self.numFaces):
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0,self.faces[id,0])
        triangle.GetPointIds().SetId(1,self.faces[id,1])
        triangle.GetPointIds().SetId(2,self.faces[id,2])
        triangles.InsertNextCell(triangle)
    self.my_vtk_polydata.SetPolys(triangles)


    curv = vtk.vtkCurvatures()
    curv.SetInputData(self.my_vtk_polydata)
    curv.SetCurvatureTypeToMean()
    curv.Update()
    pd = curv.GetOutput()
    mcurve = vtk_to_numpy((pd.GetPointData().GetArray(0)))
    
    curv.SetCurvatureTypeToGaussian()
    curv.Update()
    pd = curv.GetOutput()
    gcurve = vtk_to_numpy((pd.GetPointData().GetArray(0)))
    
    m = self.ms.current_mesh()
    density = m.vertex_custom_scalar_attribute_array("Density")
    pressure = m.vertex_custom_scalar_attribute_array("Pressure") 

    data_type = vtk.VTK_DOUBLE
    vtk_density = numpy_support.numpy_to_vtk(num_array=density, deep=True, array_type=data_type)
    vtk_density.SetName("Density")
    self.my_vtk_polydata.GetPointData().AddArray(vtk_density)

    vtk_pressure = numpy_support.numpy_to_vtk(num_array=pressure, deep=True, array_type=data_type)
    vtk_pressure.SetName("Pressure")
    self.my_vtk_polydata.GetPointData().AddArray(vtk_pressure)

    vtk_gaussian_curvature = numpy_support.numpy_to_vtk(num_array=gcurve, deep=True, array_type=data_type)
    vtk_gaussian_curvature.SetName("Gaussian Curvature")
    self.my_vtk_polydata.GetPointData().AddArray(vtk_gaussian_curvature)

    vtk_mean_curvature = numpy_support.numpy_to_vtk(num_array=mcurve, deep=True, array_type=data_type)
    vtk_mean_curvature.SetName("Mean Curvature")
    self.my_vtk_polydata.GetPointData().AddArray(vtk_mean_curvature)  
    
    return np.max(mcurve)

  
      #generate_surface_reconstruction_screened_poisson
      #compute_curvature_and_color_apss_per_vertex
      #compute_curvature_principal_directions_per_vertex#
      #compute_scalar_by_discrete_curvature_per_vertex
      #vertex_curvature_principal_dir1_matrix(self: pmeshlab.Mesh) 

  
 ################################################################# 
 # Helper functions  

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
        stepSet.append(timeStep(i,coords, i*480, (i+1)*480, "HeidiTimeStep"))

def load_scalars(filename, stepSet):
    File_data = np.loadtxt(filename , dtype=float,comments="#", skiprows=0, usecols=(5,6))
    pressure = File_data[:,0]
    density = File_data[:,1]
    #FIGURE OUT SLICING!!
    i=0
    for step in stepSet:
        #print("Step ", i, " scalar add")
        tsP = pressure[i*480:(i+1)*480]
        tsD = density[i*480:(i+1)*480]
        i=i+1
        #step.ms.print_status()
        #print("Scalars ", tsP.size, " Vertices ", step.ms.current_mesh().vertex_number())
        step.ms.current_mesh().add_vertex_custom_scalar_attribute(tsP, "Pressure")
        step.ms.current_mesh().add_vertex_custom_scalar_attribute(tsD, "Density")


######################################################################################
# MAIN SCRIPT
stepSet = []
print("Start it up")
# Text file data converted to integer data type
File_data = np.loadtxt("HeidiFieldLineXYZ2.txt", dtype=float,comments="#", skiprows=0, usecols=(2,3,4))

field_lines = np.array_split(File_data,File_data.shape[0]/5)

print("Processing ", len(field_lines), " field lines")
minb_pts=np.zeros((len(field_lines),3))
line_num=0
for line in field_lines:
    minb_pts[line_num] = line[2]
    line_num=line_num+1

find_dist_z0(minb_pts)

create_timesteps(minb_pts,540,stepSet)

load_scalars("h-clean.txt",stepSet)

i=0
maxKTimeSeries =[]
for step in stepSet:
   maxk = step.computeCurvature()
   maxz = step.get_max_z()
   maxKTimeSeries.append([i,maxk,maxz])
   step.exportSurfaceDataToVTK(step.name)
   i=i+1

# opening the csv file in 'w+' mode
file = open('meank.csv', 'w+', newline ='')
 
# writing the data into the file
with file:    
    write = csv.writer(file)
    write.writerows(maxKTimeSeries)