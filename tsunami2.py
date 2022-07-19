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
import time as time
import fem   as fem

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

def compute(theMeshFile, theResultFiles, U, V, E, dt, nIter, nSave):
  theProblem = fem.Tsunami(theMeshFile)
  theProblem.initialize(U, V, E)
  tm = time.time()
  for i in range(nIter):
    tsunamiResult(theProblem, dt)
    E = theProblem.E.reshape(theProblem.mesh.nElem, 3)
    print("Iteration = %d" % (i+1) )
    Bang = 27
    print(" == Elevations for element %d : %14.7e %14.7e %14.7e " % (Bang,*E[Bang][:]) )
    if (i + 1) % nSave == 0:
      print("Time = %f" % (time.time() - tm))   
      writeResult(theResultFiles,i + 1,E) # a revoir
      # Ne pas oublier le temps
      
  U = theProblem.U.reshape(theProblem.mesh.nElem, 3)
  V = theProblem.V.reshape(theProblem.mesh.nElem, 3)

  return [U, V, E]

# -------------------------------------------------------------------------
# TSUNAMIRESULT fait appel aux fonctions adequates pour generer les integrales de ligne et de surface
#               assemble la matrice de rigidite avec les valeurs des integrales
#               effectue une elimination de Gauss-Legendre afin d'optimiser la resolution du systeme
#               resout le systeme en tout point du maillage pour en definir les valeurs d'elevation et de vitesses
#               attend que le pas de temps soit ecoule
# @param    : theProblem        Objet representant le tsunami comme decrit dans fem.py
#             dt                Pas de temps entre les iterations
# @return   : -

def tsunamiResult(theProblem,dt):
  theProblem.FU.fill(0)
  theProblem.FV.fill(0)
  theProblem.FH.fill(0)

  addTrianglesIntegral(theProblem)
  addEdgesIntegral(theProblem)
  multiplyInverseMatrix(theProblem)

  for i in range(theProblem.size):
      theProblem.E[i] += dt * theProblem.FH[i]
      theProblem.U[i] += dt * theProblem.FU[i]
      theProblem.V[i] += dt * theProblem.FV[i]

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
      f = 2 * omega * np.sin(y/raplong)   # voir code matthieu   
      
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
          betastar = 0 #!!!!!!
          estar    = eleft + np.sqrt(h / g) * betaleft
      else:
          betastar = (betaleft + betaright)/ 2  + np.sqrt(g / h) * (eleft - eright) / 2 
          estar    = (eleft + eright) / 2       + np.sqrt(h / g) * (betaleft - betaright) / 2
          
      BH[mapEdgeLeft]  -= ((h*betastar*raplong) @ phi) * jac*wgt
      BH[mapEdgeRight] += ((h*betastar*raplong) @ phi) * jac*wgt
      
      BV[mapEdgeLeft]  -= ((ny*g*estar*raplong) @ phi) * jac*wgt
      BV[mapEdgeRight] += ((ny*g*estar*raplong) @ phi) * jac*wgt

      BU[mapEdgeLeft]  -= ((nx*g*estar*raplong) @ phi) * jac*wgt
      BU[mapEdgeRight] += ((nx*g*estar*raplong) @ phi) * jac*wgt

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
