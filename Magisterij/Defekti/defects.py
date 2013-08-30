#!/usr/bin/env python

from PyQt4.QtCore import QSize, Qt, QPointF
from PyQt4.QtGui import QApplication, QPainter, QPen, QBrush
from PyQt4.QtSvg import QSvgGenerator

import random
import math
import sys
import subprocess

IMAGE_SIZE = 400
NUM_MOLECULES = 600
MOLECULE_SIZE = 6
ARROW_SIZE = 160

def draw_molecule(painter, x, y, q, length):
  if q:
    angle = math.atan2(y - IMAGE_SIZE/2,x - IMAGE_SIZE/2) * q
  else:
    angle = math.pi / 4
  
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
  
  p = QPen(Qt.cyan)
  p.setWidth(3)
  painter.setPen(p)
  painter.drawEllipse(QPointF(IMAGE_SIZE/2, IMAGE_SIZE/2), 100, 100)
  
  
  p = QPen(Qt.blue)
  p.setWidth(5)
  painter.setPen(p)
  for i in range(12):
    angle = i * 2 * math.pi / 12
    draw_molecule(painter, IMAGE_SIZE/2 + 100 * math.cos(angle), IMAGE_SIZE/2 + 100 * math.sin(angle), q, 3 * MOLECULE_SIZE)
    
  painter.setPen(oldpen)
    
def draw_image(name, q):
  image = QSvgGenerator()
  image.setSize(QSize(IMAGE_SIZE, IMAGE_SIZE))
  image.setFileName(name + ".svg")

  painter = QPainter()
  painter.begin(image)
  painter.setBrush(Qt.white)
  painter.setPen(Qt.NoPen)
  painter.drawEllipse(QPointF(IMAGE_SIZE/2, IMAGE_SIZE/2), IMAGE_SIZE/2, IMAGE_SIZE/2)
  painter.setBrush(QBrush())
  painter.setPen(QPen())
  
  draw_defect(painter, q)
  draw_circle(painter, q)
  
  pen = QPen()
  pen.setWidth(7)
  pen.setColor(Qt.red)
  painter.setPen(pen)
  
  painter.drawLine(IMAGE_SIZE/2 - ARROW_SIZE, IMAGE_SIZE/2, IMAGE_SIZE/2 + ARROW_SIZE, IMAGE_SIZE/2)
  painter.drawLine(IMAGE_SIZE/2 + ARROW_SIZE, IMAGE_SIZE/2, IMAGE_SIZE/2 + ARROW_SIZE - 30, IMAGE_SIZE/2 + 20)
  painter.drawLine(IMAGE_SIZE/2 + ARROW_SIZE, IMAGE_SIZE/2, IMAGE_SIZE/2 + ARROW_SIZE - 30, IMAGE_SIZE/2 - 20)
  
  font = painter.font()
  font.setPixelSize(40)
  font.setBold(True)
  painter.setFont(font)
  painter.drawText(QPointF(IMAGE_SIZE/2 + ARROW_SIZE - 30, IMAGE_SIZE/2 - 30), "E")
  
  painter.end()

if __name__ == "__main__":
  app = QApplication(sys.argv)
  for i in [-1, -0.5, 0.5, 1, 1.5, 2, -2, -1.5, 0]:
    draw_image("g_defect_light_%g" % (2*i), i)
    subprocess.call(["inkscape", "--export-pdf=g_defect_light_%g.pdf" % (2*i), "g_defect_light_%g.svg" % (2*i), "--export-area-drawing"])
