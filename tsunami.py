# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Canevas de depart
# Vincent Legat
# -------------------------------------------------------------------------
# GROUPE 106
# Thomas Flamand 27721600
# Pierre Lamotte 65441500

import numpy as np
import fem as fem
#
# Fichier fem fourni par V. LEGAT pour le devoir 7
# Modifications apportées par T. FLAMAND et P. LAMOTTE
#
_gaussTri3Xsi    = np.array([1.0/6.0,2.0/3.0,1.0/6.0])
_gaussTri3Eta    = np.array([1.0/6.0,1.0/6.0,2.0/3.0])
_gaussTri3Weight = np.array([1.0/6.0,1.0/6.0,1.0/6.0])

_gaussEdg2Xsi    = np.array([-0.5773502691896257, 0.5773502691896257])
_gaussEdg2Weight = np.array([1.0,1.0])

class IntegrationRule(object):

  def __init__(self,elementType,n):
    if (elementType == "Triangle" and n == 3) :
      self.name = "Gauss with 3 points"
      self.n = n
      self.xsi    = _gaussTri3Xsi
      self.eta    = _gaussTri3Eta
      self.weight = _gaussTri3Weight
      self.phi    = [1-_gaussTri3Xsi-_gaussTri3Eta, _gaussTri3Xsi, _gaussTri3Eta]
    elif (elementType == "Edge" and n == 2) :
      self.name = "Gauss with 2 points"
      self.n = n
      self.xsi    = _gaussEdg2Xsi
      self.weight = _gaussEdg2Weight
      self.phi    = [(1 - _gaussEdg2Xsi) / 2,(1 + _gaussEdg2Xsi) / 2]
    else :
      self.name = "Unknown rule"
      self.n = 0

  def printf(self):
    print(" Integration rule : %s " % self.name)
    print(" Number of nodes = %d " % self.n)
    print(" xsi     = ",self.xsi)
    print(" eta     = ",self.eta)
    print(" weights = ",self.weight)

# -------------------------------------------------------------------------

class Mesh(object):

  def __init__(self,filename):
    with open(filename,"r") as f :
      self.nNode = int(f.readline().split()[3])
      self.xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(self.nNode)))
      self.nElem = int(f.readline().split()[3])
      self.elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(self.nElem)))
      self.X     = self.xyz[:,0]
      self.Y     = self.xyz[:,1]
      self.H     = self.xyz[:,2]


  def printf(self):
    print("Number of nodes %d" % self.nNode)
    for i in range(self.nNode):
      print("%6d : %14.7e %14.7e" % (i,*self.xyz[i,:]))
    print("Number of triangles %d" % self.nElem)
    for i in range(self.nElem):
      print("%6d : %6d %6d %6d" % (i,*self.elem[i,:]))

  def write(self,filename,iter,E):
    fileName = fileBaseName % iter
    nElem = E.shape[0]
    with open(filename,"w") as f :
      f.write("Number of elements %d\n" % nElem)
      for i in range(nElem):
        f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*E[i,:]))
      print(" === iteration %6d : writing %s ===" % (iter,fileName))

# -------------------------------------------------------------------------

class Edges(object):

  def __init__(self,mesh):
    self.mesh = mesh
    self.nEdges = mesh.nElem * 3
    self.nBoundary = 0
    self.edges = [[0 for i in range(4)] for i in range(self.nEdges)]
    for i in range (mesh.nElem) :
      for j in range(3) :
        id = i*3 + j
        self.edges[id][0] = mesh.elem[i][j]
        self.edges[id][1] = mesh.elem[i][(j+1)%3]
        self.edges[id][2] = i
        self.edges[id][3] = -1
    self.edges.sort(key = lambda item : -(min(item[0:2])*self.nEdges)-max(item[0:2]))
    index = 0
    for i in range(self.nEdges) :
      if (self.edges[i][0:2] != self.edges[i-1][1::-1]) :
         self.edges[index] = self.edges[i]
         index += 1
      else :
         self.edges[index-1][3] = self.edges[i][2]
    del self.edges[index:]
    self.edges.sort(key = lambda item : item[3])
    self.nBoundary = 2*index - self.nEdges
    self.nEdges = index

  def printf(self):
    print("Number of edges %d" % self.nEdges)
    print("Number of boundary edges %d" % self.nBoundary)
    for i in range(self.nEdges):
      print("%6d : %4d %4d : %4d %4d" % (i,*self.edges[i]))

# -------------------------------------------------------------------------

class Tsunami(object):

  def __init__(self, filename):
    self.mesh = Mesh(filename)
    self.edges = Edges(self.mesh)
    self.rule = IntegrationRule("Triangle", 3)
    self.ruleEdge = IntegrationRule("Edge",2)

    size = 3 * self.mesh.nElem
    self.size = size
    self.E = np.zeros(size)
    self.U = np.zeros(size)
    self.V = np.zeros(size)
    self.FU = np.zeros(size)
    self.FV = np.zeros(size)
    self.FH = np.zeros(size)

    self.mapEdgeLeft = np.zeros((self.edges.nEdges,2),dtype=np.int)
    self.mapEdgeRight = np.zeros((self.edges.nEdges,2),dtype=np.int)
    for iEdge in range(self.edges.nBoundary,self.edges.nEdges):
      myEdge = self.edges.edges[iEdge]
      elementLeft  = myEdge[2]
      elementRight = myEdge[3]
      nodesLeft    = self.mesh.elem[elementLeft]
      nodesRight   = self.mesh.elem[elementRight]
      self.mapEdgeLeft[iEdge,:]  = [3*elementLeft  + np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
      self.mapEdgeRight[iEdge,:] = [3*elementRight + np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]

    self.mapTriangle = np.zeros((self.mesh.nElem,3),dtype=np.int)
    for iElem in range(self.mesh.nElem) :
      self.mapTriangle[iElem,:] = [3*iElem+j for j in range(3)]

# -------------------------------------------------------------------------

  def initialize(self,Ui,Vi,Ei):

    self.E = Ei.reshape(len(Ei)*3)
    self.U = Ui.reshape(len(Ui)*3)
    self.V = Vi.reshape(len(Vi)*3)

#
# Fin Fichier fem 
#
    
g = 9.81
R = 6371220
gamma = 10**(-7)
omega = 2 * np.pi / 86400

def readMesh(fileName) :
  with open(fileName,"r") as f :
    nNode = int(f.readline().split()[3])
    xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(nNode)))
    nElem = int(f.readline().split()[3])
    elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(nElem)))
  X = xyz[:,0]
  Y = xyz[:,1]
  H = xyz[:,2]
  return [nNode,X,Y,H,nElem,elem]

# -------------------------------------------------------------------------

def readResult(fileBaseName,iter,nElem) :
  fileName = fileBaseName % iter
  with open(fileName,"r") as f :
    nSize = int(f.readline().split()[3])
    if (nElem != nSize) :
      print(" ==== Error : incoherent sizes : %d != %d" % (nElem,nSize))
    E = np.array(list(list(float(w) for w in f.readline().split()[2:5]) for i in range(nElem)))
    print(" === iteration %6d : reading %s ===" % (iter,fileName))
  # Toujours fermer un ficher ouvert
  f.close()
  return E

# -------------------------------------------------------------------------

def writeResult(fileBaseName,iter,E) :
  fileName = fileBaseName % iter
  nElem = E.shape[0]
  with open(fileName,"w") as f :
    f.write("Number of elements %d\n" % nElem)
    for i in range(nElem):
      f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*E[i,:]))
    print(" === iteration %6d : writing %s ===" % (iter,fileName))

# -------------------------------------------------------------------------

def initialConditionOkada(x,y) :
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

# -------------------------------------------------------------------------
# Debut de la section de code a realiser
# Rappel des fonctions de formes :
#       2D -> phi = [(1 - xsi)/2,(1 + xsi)/2]
#       3D -> phi = [1 - xsi - eta,xsi,eta]

# -------------------------------------------------------------------------
# COMPUTE lit le ficher contenant le maillage et fait appel aux fonctions de resolution
#         ecrit l'elevation de l'ocean en tout point dans des fichiers de sauvegarde
#         retourne les matrices de vitesses et d'elevation sur tout le maillage
# @param    : theMeshFile       Fichier contenant le maillage a analyser
#             theResultFiles    Fichiers d'ecriture des resultats d'elevation
#             U                 Vecteur vitesse selon x de la colonne d'eau en tout point
#             V                 Vecteur vitesse selon y de la colonne d'eau en tout point
#             E                 Matrice des conditions initiales pour chaque triangle
#             dt                Pas de temps entre les iterations
#             nIter             Nombre d'iterations
#             nSave             Nombre de fichiers de sauvegarde
# @return   : U                 Vecteur vitesse selon x de la colonne d'eau en tout point
#             V                 Vecteur vitesse selon y de la colonne d'eau en tout point
#             E                 Matrice des conditions initiales pour chaque triangle

def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):
    theProblem = fem.Tsunami(theMeshFile)
    theProblem.initialize(U, V, E)
    for i in range(nIter): 
      theProblem.FU.fill(0)
      theProblem.FV.fill(0)
      theProblem.FH.fill(0)
       
      addTrianglesIntegral(theProblem)
      addEdgesIntegral(theProblem)
      multiplyInverseMatrix(theProblem)
       
      theProblem.E += dt * theProblem.FH
      theProblem.U += dt * theProblem.FU
      theProblem.V += dt * theProblem.FV

        
      if np.mod(i+1,nSave) == 0:
          writeResult(theResultFiles,i+1,E.reshape(theProblem.mesh.nElem, 3))
             
    return [theProblem.U.reshape(theProblem.mesh.nElem, 3),theProblem.V.reshape(theProblem.mesh.nElem, 3),theProblem.E.reshape(theProblem.mesh.nElem, 3)]   


# -------------------------------------------------------------------------
# ADDTRIANGLESINTEGRAL genere la valeur des integrales de surface
# @param    : theProblem        Objet representant le tsunami comme decrit dans fem.py
# @return   :

def addTrianglesIntegral(theProblem):
  theMesh = theProblem.mesh
  theRule = theProblem.rule
  U  = theProblem.U
  V  = theProblem.V
  E  = theProblem.E
  X = theMesh.X
  Y = theMesh.Y
  H = theMesh.H  
  BU = theProblem.FU
  BV = theProblem.FV
  BH = theProblem.FH
  phi = theRule.phi
  wgt = theRule.weight

  for iElem in range(theMesh.nElem):
      [dphidx, dphidy, jac] = computeShapeTriangle(theMesh, iElem)
      mapTriangle = theProblem.mapTriangle[iElem]
      mapCoord = theMesh.elem[iElem]
      
      u = phi @ U[mapTriangle]
      v = phi @ V[mapTriangle]
      e = phi @ E[mapTriangle]        
      x = phi @ X[mapCoord]
      y = phi @ Y[mapCoord]
      h = phi @ H[mapCoord]
      
      raplong = (4*R**2 + x**2 + y**2)/(4*R**2) # raplong = X/X* = Y/Y*
      f = 2 * omega * (4*R*R - x*x - y*y) / (4*R*R + x*x + y*y)           
      
      BH[mapTriangle] += sum( (np.outer(h*u*raplong,dphidx) + np.outer(h*v*raplong,dphidy))* jac*wgt)  + (((h*x*u)/(R*R)  + (h*y*v)/(R*R)) @ phi )* jac*wgt
      BV[mapTriangle] +=  (((-f*u) - (gamma*v) + ((g*y*e) / (2*R*R))) @ phi)*jac*wgt + sum(np.outer(g*e*raplong,dphidy)* jac*wgt)
      BU[mapTriangle] +=  (((f*v) - (gamma*u) + ((g*x*e) / (2*R*R))) @ phi)*jac*wgt + sum(np.outer(g*e*raplong,dphidx)* jac*wgt)

# -------------------------------------------------------------------------
# ADDEDGESINTEGRAL genere la valeur des integrales de ligne
# @param    : theProblem        Objet representant le tsunami comme decrit dans fem.py
# @return   : -

def addEdgesIntegral(theProblem):
  theMesh  = theProblem.mesh
  theRule  = theProblem.ruleEdge
  theEdges = theProblem.edges
  U   = theProblem.U
  V   = theProblem.V
  E   = theProblem.E
  H = theMesh.H
  X = theMesh.X
  Y = theMesh.Y  
  BU  = theProblem.FU
  BV  = theProblem.FV
  BH  = theProblem.FH
  phi = theRule.phi
  wgt = theRule.weight

  for iEdge in range(theEdges.nEdges):
      [phiphi, nx, ny, jac] = computeShapeEdge(theEdges, iEdge)
      mapEdgeLeft  = theProblem.mapEdgeLeft[iEdge]
      mapEdgeRight = theProblem.mapEdgeRight[iEdge]
      EdgeCoord = theEdges.edges[iEdge][0:2]
      
      uleft  = phi @ U[mapEdgeLeft]
      vleft  = phi @ V[mapEdgeLeft]
      uright = phi @ U[mapEdgeRight]
      vright = phi @ V[mapEdgeRight]
      eleft  = phi @ E[mapEdgeLeft]
      eright = phi @ E[mapEdgeRight]
      betaleft  = uleft * nx  + vleft * ny
      betaright = uright * nx + vright * ny
      
      h = phi @ H[EdgeCoord]
      x = phi @ X[EdgeCoord]
      y = phi @ Y[EdgeCoord]
          
      raplong = (4*R**2 + x**2 + y**2)/(4*R**2) # raplong = X/X* = Y/Y*          
          
      if theEdges.edges[iEdge][3] == -1:
          betastar = 0 
          estar    = eleft + np.sqrt(h / g) * betaleft
      else:
          betastar = (betaleft + betaright)/ 2  + np.sqrt(g / h) * (eleft - eright) / 2 
          estar    = (eleft + eright) / 2       + np.sqrt(h / g) * (betaleft - betaright) / 2
          
      BH[mapEdgeLeft]  -= ((h*betastar*raplong) @ phi) * jac*wgt/2
      BH[mapEdgeRight] += ((h*betastar*raplong) @ phi) * jac*wgt/2
      
      BV[mapEdgeLeft]  -= ((ny*g*estar*raplong) @ phi) * jac*wgt/2
      BV[mapEdgeRight] += ((ny*g*estar*raplong) @ phi) * jac*wgt/2

      BU[mapEdgeLeft]  -= ((nx*g*estar*raplong) @ phi) * jac*wgt/2
      BU[mapEdgeRight] += ((nx*g*estar*raplong) @ phi) * jac*wgt/2

# -------------------------------------------------------------------------

def multiplyInverseMatrix(theProblem):
  theMesh  = theProblem.mesh
  BU = theProblem.FU
  BV = theProblem.FV
  BH = theProblem.FH
  Ainverse = np.array([[18.0, -6.0, -6.0], [-6.0, 18.0, -6.0], [-6.0, -6.0, 18.0]])
  for iElem in range(theMesh.nElem) :
     map   = theProblem.mapTriangle[iElem]
     nodes = theMesh.elem[iElem]
     x   = theMesh.X[nodes]
     y   = theMesh.Y[nodes]
     jac = abs((x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]))
     BU[map] = Ainverse @ BU[map] / jac
     BV[map] = Ainverse @ BV[map] / jac
     BH[map] = Ainverse @ BH[map] / jac

# -------------------------------------------------------------------------

def computeShapeTriangle(theMesh,theElement) :
  dphidxsi = np.array([ 1.0, 0.0,-1.0])
  dphideta = np.array([ 0.0, 1.0,-1.0])

  nodes = theMesh.elem[theElement]
  x = theMesh.X[nodes]
  y = theMesh.Y[nodes]
  dxdxsi = x @ dphidxsi
  dxdeta = x @ dphideta
  dydxsi = y @ dphidxsi
  dydeta = y @ dphideta
  jac    = abs(dxdxsi * dydeta - dxdeta*  dydxsi)
  dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jac
  dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jac
  return [dphidx, dphidy, jac]

# -------------------------------------------------------------------------

def computeShapeEdge(theEdges, iEdge):
  phiphi   = np.array([[ 2.0, 1.0],[ 1.0, 2.0]]) / 3.0

  nodes = theEdges.edges[iEdge][0:2]
  x   = theEdges.mesh.X[nodes]
  y   = theEdges.mesh.Y[nodes]
  dx  = x[1] - x[0]
  dy  = y[1] - y[0]
  jac = np.sqrt(dx * dx + dy * dy)
  nx  =  dy / jac
  ny  = -dx / jac
  return [phiphi, nx, ny, jac]

