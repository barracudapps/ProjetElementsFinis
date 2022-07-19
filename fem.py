#
# Fichier fem fourni par V. LEGAT pour le devoir 7
# Modifications apportées par T. FLAMAND et P. LAMOTTE
#
import numpy as np
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
