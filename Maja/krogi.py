from PyQt4.QtGui import QImage, QPainter, QApplication, QPen
from PyQt4.QtCore import Qt, QPointF

import sys
from math import sqrt

WIDTH = 1600
HEIGHT = 1200

SOURCE_ONE = 600
SOURCE_TWO = 1000
SOURCE_DIF = 400.0

POS_Y = 600
RADIJ = 70.0

POINT_SIZE = 6.0


def krogi(n):
	image = QImage(1600, 1200, QImage.Format_RGB32)
	image.fill(Qt.white)
	painter = QPainter(image)
	
	pen = QPen()
	pen.setColor(Qt.blue)
	pen.setWidth(3)
	painter.setPen(pen)
	
	for i in range(n):
		r = RADIJ * 2 * (i+1)
		painter.drawEllipse(SOURCE_ONE - r, POS_Y-r, 2*r, 2*r)
		
	for i in range(n):
		r = RADIJ * 2 * (i+1)
		painter.drawEllipse(SOURCE_TWO - r, POS_Y-r, 2*r, 2*r)
		
	
	pen = QPen()
	pen.setColor(Qt.darkGray)
	pen.setWidth(2)
	painter.setPen(pen)
	
	for i in range(n):
		r = RADIJ * (2*i+1)
		painter.drawEllipse(SOURCE_ONE - r, POS_Y-r, 2*r, 2*r)
		
	for i in range(n):
		r = RADIJ * (2*i+1)
		painter.drawEllipse(SOURCE_TWO - r, POS_Y-r, 2*r, 2*r)
		
	# Pikice
	painter.setPen(Qt.red)
	painter.setBrush(Qt.red)
	SourceOne = QPointF(SOURCE_ONE, POS_Y)
	for i in range(n):
		for j in range(n):
			# c^2 = a^2 + b^2 - 2ab cos gama
			# => cos gama = (a^2+b^2-c^2) / 2ab
			a = RADIJ * 2 * (j+1)
			b = RADIJ * 2 * (i+1)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
			a = RADIJ * (2*j+1)
			b = RADIJ * (2*i+1)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
	painter.setPen(Qt.darkGreen)
	painter.setBrush(Qt.darkGreen)
	SourceOne = QPointF(SOURCE_ONE, POS_Y)
	for i in range(n):
		for j in range(n):
			# c^2 = a^2 + b^2 - 2ab cos gama
			# => cos gama = (a^2+b^2-c^2) / 2ab
			a = RADIJ * (2*j+2)
			b = RADIJ * (2*i+1)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
			a = RADIJ * (2*j+1)
			b = RADIJ * (2*i+2)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
	painter.setPen(Qt.black)
	painter.setBrush(Qt.black)
	
	painter.drawEllipse(QPointF(SOURCE_ONE, POS_Y), 8, 8)
	painter.drawEllipse(QPointF(SOURCE_TWO, POS_Y), 8, 8)
		
	painter.end()
	image.save("krogi.png")
	
def svetla(n):
	image = QImage(1600, 1200, QImage.Format_RGB32)
	image.fill(Qt.white)
	painter = QPainter(image)
	
	POINT_SIZE = 4.0
	
	pen = QPen()
	pen.setColor(Qt.blue)
	pen.setWidth(1)
	painter.setPen(pen)
	
	for i in range(n):
		r = RADIJ * 2 * (i+1)
		painter.drawEllipse(SOURCE_ONE - r, POS_Y-r, 2*r, 2*r)
		
	for i in range(n):
		r = RADIJ * 2 * (i+1)
		painter.drawEllipse(SOURCE_TWO - r, POS_Y-r, 2*r, 2*r)
		
	
	pen = QPen()
	pen.setColor(Qt.darkGray)
	pen.setWidth(1)
	painter.setPen(pen)
	
	for i in range(n):
		r = RADIJ * (2*i+1)
		painter.drawEllipse(SOURCE_ONE - r, POS_Y-r, 2*r, 2*r)
		
	for i in range(n):
		r = RADIJ * (2*i+1)
		painter.drawEllipse(SOURCE_TWO - r, POS_Y-r, 2*r, 2*r)
		
	# Pikice
	painter.setPen(Qt.red)
	painter.setBrush(Qt.red)
	SourceOne = QPointF(SOURCE_ONE, POS_Y)
	for i in range(n):
		for j in range(n):
			# c^2 = a^2 + b^2 - 2ab cos gama
			# => cos gama = (a^2+b^2-c^2) / 2ab
			a = RADIJ * 2 * (j+1)
			b = RADIJ * 2 * (i+1)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
			a = RADIJ * (2*j+1)
			b = RADIJ * (2*i+1)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
	painter.setPen(Qt.darkGreen)
	painter.setBrush(Qt.darkGreen)
	SourceOne = QPointF(SOURCE_ONE, POS_Y)
	for i in range(n):
		for j in range(n):
			# c^2 = a^2 + b^2 - 2ab cos gama
			# => cos gama = (a^2+b^2-c^2) / 2ab
			a = RADIJ * (2*j+2)
			b = RADIJ * (2*i+1)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
			a = RADIJ * (2*j+1)
			b = RADIJ * (2*i+2)
			c = SOURCE_DIF
						
			cA = (c*c+b*b-a*a) / 2 / b / c
			if abs(cA) < 1:
				print(', '.join([str(i) for i in [a, b, c, cA]]))
				sA = sqrt(1-cA*cA)
				p = SourceOne + QPointF(b * cA, b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				p = SourceOne + QPointF(b * cA, -b * sA)
				painter.drawEllipse(p, POINT_SIZE, POINT_SIZE)
				
	painter.setPen(Qt.black)
	painter.setBrush(Qt.black)
	
	painter.drawEllipse(QPointF(SOURCE_ONE, POS_Y), 8, 8)
	painter.drawEllipse(QPointF(SOURCE_TWO, POS_Y), 8, 8)
		
	painter.end()
	image.save("svetla.png")

if __name__ == "__main__":
	app = QApplication(sys.argv)
	krogi(8)
	svetla(8)