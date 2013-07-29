from PyQt4.QtCore import QSize, Qt, QPointF
from PyQt4.QtGui import QApplication, QPainter
from PyQt4.QtSvg import QSvgGenerator

import random
import math
import sys
import subprocess

IMAGE_SIZE = 400
NUM_MOLECULES = 600
MOLECULE_SIZE = 6

def draw_molecule(painter, x, y, q, length):
  angle = math.atan2(y - IMAGE_SIZE/2,x - IMAGE_SIZE/2) * q
  
  dx = length * math.cos(angle)
  dy = length * math.sin(angle)
  
  painter.drawLine(x-dx,y-dy,x+dx,y+dy)

def radius(x, y):
  return (x - IMAGE_SIZE/2) * (x - IMAGE_SIZE/2) + (y - IMAGE_SIZE/2) * (y - IMAGE_SIZE/2) 

def draw_defect(painter, q):
  for i in range(NUM_MOLECULES):
    x = random.random() * IMAGE_SIZE
    y = random.random() * IMAGE_SIZE
    if radius(x, y) < IMAGE_SIZE * IMAGE_SIZE / 4:
      draw_molecule(painter, x, y, q, MOLECULE_SIZE)
    
      
def draw_circle(painter, q):
  oldpen = painter.pen()
  painter.setPen(Qt.cyan)
  painter.drawEllipse(QPointF(IMAGE_SIZE/2, IMAGE_SIZE/2), 100, 100)
  
  painter.setPen(Qt.blue)
  for i in range(12):
    angle = i * 2 * math.pi / 12
    draw_molecule(painter, IMAGE_SIZE/2 + 100 * math.cos(angle), IMAGE_SIZE/2 + 100 * math.sin(angle), q, 3 * MOLECULE_SIZE)
    
def draw_image(name, q):
  image = QSvgGenerator()
  image.setSize(QSize(IMAGE_SIZE, IMAGE_SIZE))
  image.setFileName(name + ".svg")

  painter = QPainter()
  painter.begin(image)
  draw_defect(painter, q)
  draw_circle(painter, q)
  painter.end()

if __name__ == "__main__":
  app = QApplication(sys.argv)
  for i in [-1, -0.5, 0.5, 1]:
    draw_image("g_defect_%g" % i, i)
    subprocess.call(["inkscape", "--export-pdf=g_defect_%g.pdf" % i, "g_defect_%g.svg" % i, "--export-area-drawing"])
