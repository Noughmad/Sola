from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtSvg import *

import sys
import math
import random

LENGTH = 25
JITTER_MAX = math.pi / 100


if __name__ == "__main__":
  app = QApplication(sys.argv)
  
  generator = QSvgGenerator()
  generator.setFileName("radial-cross.svg")
  generator.setSize(QSize(200, 200))
  generator.setViewBox(QRect(0, 0, 200, 200))
  
  painter = QPainter()
  painter.begin(generator)
  
  pen = QPen()
  pen.setWidth(3)
  painter.setPen(pen)
  painter.drawEllipse(QPoint(100, 100), 100, 100)
  pen.setWidth(8)
  painter.setPen(pen)
  painter.drawEllipse(QPoint(100, 100), 2, 2)
  
  painter.setPen(QPen())
  
  for i in range(100):
    startAngle = 2 * math.pi * i / 100
    endAngle = startAngle + JITTER_MAX * (2 * random.random() - 1)
    startPos = (1 - math.pow(random.random(),2)) * (100 - LENGTH)
    
    painter.drawLine(
      100 + startPos * math.cos(startAngle),
      100 + startPos * math.sin(startAngle),
      100 + (startPos + LENGTH) * math.cos(endAngle),
      100 + (startPos + LENGTH) * math.sin(endAngle)
    )
  
  painter.end()