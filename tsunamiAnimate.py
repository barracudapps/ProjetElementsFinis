#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Visualisation des résultats : OpenGL old-fashion
# Une implémentation avec des shaders sera présentée en S11
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *

import numpy as np
import tsunami
import sys
import time

# -------------------------------------------------------------------------

def draw():  
  global E,theFlagBathymetry,theMouse,theRatio,translate,sky,timeZero,follow,lastFrame

  glClearColor( sky[0], sky[1], sky[2], 0.0 );
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(theRatio*65.0,theRatio,theRatio,100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if follow:
    dt = 12 # dt = 12 optimal pour les résultats fournis par le professeur
            # A modifier selon le pas de temps des fichiers
    t = time.time() - timeZero
    frame = t * 2
    if (dt * frame < 136.0):
      gluLookAt(-1.0,6.0 -(0.001*np.exp(0.063*frame*dt)),1.0,0.0,15.0,5.0 -(frame*dt*0.038),0.0,0.0,1.0);
      lastFrame = t
    else:
      gluLookAt(-1.0,6.0 -(0.001*np.exp(0.063*2*lastFrame*dt)),1.0,0.0,15.0,5.0 -(2*lastFrame*dt*0.038),0.0,0.0,1.0);
  else:
    gluLookAt(0.0,1.0,0.0,0.0,20.0,0.0,0.0,0.0,1.0);  
  glTranslatef(0.0,14.0,translate+0.0);
  glRotatef(0.3*theMouse,0.0,0.0,1.0);
  
  quadratic = gluNewQuadric();         
  gluQuadricNormals(quadratic, GLU_SMOOTH); 
  glColor3f(0.4,0.7,0.5);
  gluSphere(quadratic,5.95,400,200);
 
  n = 9*nElem
  index  = np.ravel(elem)
  colors = np.zeros(n) 
  coord  = np.zeros(n)

  if (theFlagBathymetry == False) :
    value = 10* np.ravel(E)
  else :
    value = H[index]/BathMax
    
  value = np.clip(value,0,1)        
  colors[0:n:3] = 3.5*(value)*(value)
  colors[1:n:3] = (1-value)*(value)*3.5
  colors[2:n:3] = (1-value)*(1-value)
  
  x = X[index]
  y = Y[index]
  factor = (4*R*R + x*x + y*y)*(R/6)
  coord[0:n:3] = 4*R*R * x / factor
  coord[1:n:3] = 4*R*R * y / factor
  coord[2:n:3] = (4*R*R - x*x - y*y)*R / factor   
 
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_FLOAT,0,coord);
  glNormalPointer(GL_FLOAT,0,coord);
  glColorPointer(3,GL_FLOAT,0,colors);
  glDrawArrays(GL_TRIANGLES,0,nElem*3);
  glDisableClientState(GL_NORMAL_ARRAY);    
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);      
    
  coord *= 1.001;

  glColor3f(0.0, 0.0, 0.0);
  glEnableClientState(GL_VERTEX_ARRAY);    
  glVertexPointer(3,GL_FLOAT,0,coord);
  for i in range(nElem):
    glDrawArrays(GL_LINE_LOOP,3*i,3)
  glDisableClientState(GL_VERTEX_ARRAY) 
  glutSwapBuffers() 
 
# -------------------------------------------------------------------------
    
def reshape(width, height):   
  global theRatio
  
  glViewport(0,0,width,height)
  theRatio = width / float(height)
  glutPostRedisplay()

# -------------------------------------------------------------------------
 
def keyboard(key,x,y):
  global theFlagBathymetry,pause,restart,zoom,theRatio,sky,follow
  
  key = key.decode()
  if ord(key) == 27: # Escape
    sys.exit(0)
  elif key == 'b':
    theFlagBathymetry = True
  elif key == 'c':
    if (sky == [0.9, 0.9, 0.8]):
      sky = [1.0, 1.0, 1.0]
    elif (sky == [1.0, 1.0, 1.0]):
      sky = [0.2, 0.2, 0.2]
    elif (sky == [0.2, 0.2, 0.2]):
      sky = [0.9, 0.9, 0.8]
  elif key == 'e':
    theFlagBathymetry = False
  elif key == 'f':
    follow = not follow
  elif key == 'p':
    pause = not pause
  elif key == 'r':
    restart = not restart
  elif key == 'z':
    if zoom:
      theRatio = theRatio * 2
      zoom = False
    else:
      theRatio = theRatio / 2
      zoom = True
  else:
    return
  glutPostRedisplay()

# -------------------------------------------------------------------------

def special(symbol,x,y):
  global theMouse,translate
  
  if symbol == GLUT_KEY_LEFT :        
    theMouse -= 5
  elif symbol == GLUT_KEY_RIGHT :
    theMouse += 5
  elif symbol == GLUT_KEY_UP :
    translate += 0.5
  elif symbol == GLUT_KEY_DOWN :
    translate -= 0.5
  else:
    return
  glutPostRedisplay()
  
# -------------------------------------------------------------------------

def idle():
  global iter,delta,E,theResultFiles,restart,pause,timeZero

  if restart:
    iter = 0
    timeZero = time.time()-0.35
    restart = False
    pause = False
  if(not pause):
    iter += delta 
    try :
      E = tsunami.readResult(theResultFiles,iter,nElem)
      glutPostRedisplay()
    except FileNotFoundError: 
      pass

# -------------------------------------------------------------------------
  
iter = 0; delta = 25;
R = 6371220;
BathMax = 9368;
theMeshFile = "PacificTriangleFine.txt"
theResultFiles = "results/eta-%06d.txt"
theFlagBathymetry = False
translate = 0
theMouse = 389
theRatio = 1.0
restart = False
pause = False
zoom = False
playMusic = True
follow  = True
lastFrame = 0.0

sky = [0.9, 0.9, 0.8]

glutInit(sys.argv)
timeZero = time.time()
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
glutInitWindowPosition(100, 100)
glutInitWindowSize(500, 500)
glutCreateWindow(b"MECA1120 : the 2019 project :-)")

matSpecular   = [1.0,1.0,1.0,0.0]
matShininess  = [50.0]
lightPosition = [8.0,8.0,8.0,0.0]
lightRadiance = [1.0,1.0,1.0,1.0]
glMaterialfv(GL_FRONT,GL_SPECULAR, matSpecular)
glMaterialfv(GL_FRONT,GL_SHININESS,matShininess)
glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE)
glLightfv(GL_LIGHT0,GL_POSITION,lightPosition)
glLightfv(GL_LIGHT0,GL_DIFFUSE, lightRadiance)
glLightfv(GL_LIGHT0,GL_SPECULAR,lightRadiance)
glEnable(GL_LIGHTING)
glEnable(GL_LIGHT0)
glDepthFunc(GL_LEQUAL)
glEnable(GL_DEPTH_TEST)
glEnable(GL_COLOR_MATERIAL)
glEnable(GL_NORMALIZE)	

glutDisplayFunc(draw)
glutKeyboardFunc(keyboard)
glutSpecialFunc(special)
glutReshapeFunc(reshape)
glutIdleFunc(idle)

glClearColor( 0.9, 0.9, 0.8, 0.0 );
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
gluPerspective(65.0,1.0,1.0,100.0);

if "-info" in sys.argv:
  print("GL_RENDERER   = ",glGetString(GL_RENDERER).decode())
  print("GL_VERSION    = ",glGetString(GL_VERSION).decode())
  print("GL_VENDOR     = ",glGetString(GL_VENDOR).decode())

print('======================================')  
print(' b          : show bathymetry ')
print(' c          : change background color ')
print(' e          : show elevation (by default) ')
print(' f          : follow the wave ')
print(' p          : pause/start animation ')
print(' r          : restart animation ')
print(' z          : zoom in/out the Earth ')
print(' UP/DOWN    : translate the Earth ')
print(' LEFT/RIGHT : rotate the Earth ')
print(' ESC        : exit ')
print('======================================')
 
[nNode,X,Y,H,nElem,elem] = tsunami.readMesh(theMeshFile)
try :
  E = tsunami.readResult(theResultFiles,0,nElem)
except FileNotFoundError:
  E = np.zeros([nElem,3])

glutMainLoop()

# -------------------------------------------------------------------------
