import numpy as np
import fem as fem
import time as time

def initialConditionOkada(x,y) :
  R = 6371220;
  x3d = 4*R*R*x / (4*R*R + x*x + y*y);
  y3d = 4*R*R*y / (4*R*R + x*x + y*y);
  z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  lat = np.arcsin(z3d/R)*180/np.pi;
  lon = np.arctan2(y3d,x3d)*180/np.pi;
  lonMin = 142;
  lonMax = 143.75;
  latMin = 35.9;
  latMax = 39.5;
  olon = (lonMin+lonMax)/2;
  olat = (latMin+latMax)/2;
  angle = -12.95*np.pi/180; 
  lon2 = olon + (lon-olon)*np.cos(angle) + (lat-olat)*np.sin(angle);
  lat2 = olat - (lon-olon)*np.sin(angle) + (lat-olat)*np.cos(angle);
  return np.all([lon2 <= lonMax,lon2 >= lonMin,lat2 >= latMin,lat2 <= latMax],axis=0).astype(int)
 

theMeshFile = "PacificTiny.txt"
theProblem = fem.Tsunami(theMeshFile)
theMesh = theProblem.mesh
theEdges = theProblem.edges

#nElem = theMesh.nElem
#elem = theMesh.elem
#X = theMesh.X
#Y = theMesh.Y
#x = np.zeros([nElem,3])
#y = np.zeros([nElem,3])
#for iElem in range(nElem):
#  nodes  = elem[iElem]
#  x[iElem][:] = X[nodes]
#  y[iElem][:] = Y[nodes]
#E = initialConditionOkada(x,y)
#U = np.zeros([nElem,3])
#V = np.zeros([nElem,3])
#
#
#theProblem.initialize(U,V,E)
#print(theProblem.U)
#print(theProblem.E)

#print(theMesh.X)
#print(theMesh.nElem)
#print(theEdges.nBoundary)
#print(theEdges.edges[0][3])
#U = theProblem.U
#M = np.array([[0,1,2],[3,4,5],[6,7,8],[10,11,12]])
#print(len(M))
#print(M.reshape(len(M)*3))
for iElem in range(10):
    mapTriangle = theProblem.mapTriangle[iElem]
    print(theMesh.H[iElem])
    for i in range(2):
          
        print(theMesh.H[theEdges.edges[iElem][i]])
