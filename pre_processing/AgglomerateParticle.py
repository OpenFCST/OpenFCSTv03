import smesh, geompy, SMESH
import numpy
#import SALOMEDS
import math

class AgglomerateParticle:
    
    def __init__(self, position, radius, filmThickness, filmOrient, percentageCoverage):
        self.position = position
        self.radius = radius
        self.filmThickness =filmThickness
        self.filmOrient = filmOrient
        self.percentageCoverage = percentageCoverage
        
    def generateCarbonPlatinumGeometry(self):
        self.__carbonPlatinumSphere = geompy.MakeSpherePntR(geompy.MakeVertex(self.position[0],self.position[1], self.position[2] ), self.radius)
    
    def getCarbonParticle(self):
        return self.__carbonPlatinumSphere
    
    def getIonomerPipe(self):
        return self.__ionomerPipe
    
    def generateFilmGeometry(self):
        if self.percentageCoverage < 0:
            print "Error: Cannot compute negative coverage percentages!!!"
            raise Exception
        elif self.percentageCoverage <100:
            theta=(math.pi)*float(self.percentageCoverage/100)
            zetha = theta/2 +self.filmOrient[0] -(math.pi)/2
            gamma = self.filmOrient[0] - theta/2
            print theta ,zetha, gamma
            arc1Pt1 = [-math.sin(zetha)*self.radius + self.position[0], math.cos(zetha)*self.radius + self.position[1], self.position[2]] 
            arc1Pt2 = [math.sin(gamma)*self.radius  + self.position[0], math.cos(gamma)*self.radius  + self.position[1], self.position[2] ]
            
            arc2Pt1 = [arc1Pt1[0] -math.sin(zetha)*self.filmThickness, arc1Pt1[1] +math.cos(zetha)*self.filmThickness, arc1Pt1[2] ]
            arc2Pt2 = [arc1Pt2[0] +math.sin(gamma)*self.filmThickness, arc1Pt2[1] +math.cos(gamma)*self.filmThickness, arc1Pt2[2] ]
            
            circleCentre = geompy.MakeVertex(self.position[0],self.position[1], self.position[2] )
            
            arc1V1 =  geompy.MakeVertex(arc1Pt1[0],arc1Pt1[1],arc1Pt1[2] )
            arc1V2 =  geompy.MakeVertex(arc1Pt2[0],arc1Pt2[1],arc1Pt2[2] )
            
            arc2V1 =  geompy.MakeVertex(arc2Pt1[0],arc2Pt1[1],arc2Pt1[2] )
            arc2V2 =  geompy.MakeVertex(arc2Pt2[0],arc2Pt2[1],arc2Pt2[2] )
            
            print arc1Pt1, arc1Pt2, arc2Pt1, arc2Pt2
            arc1 = geompy.MakeArcCenter(circleCentre, arc1V1, arc1V2)
            arc2 = geompy.MakeArcCenter(circleCentre, arc2V1, arc2V2) 
            
            line1 = geompy.MakeLineTwoPnt(arc1V1, arc2V1)
            line2 = geompy.MakeLineTwoPnt(arc1V2, arc2V2)
            
            face1 = geompy.MakeFaceWires([arc1, arc2, line1, line2], 1)
            
            circleVector = geompy.MakeVectorDXDYDZ(self.filmOrient[0], self.filmOrient[1], self.filmOrient[2])
            circle1 = geompy.MakeCircle(circleCentre, circleVector , self.radius)
            
            self.__ionomerPipe = geompy.MakePipe(face1, circle1)
        elif self.percentageCoverage == 100:
            circleCentre = geompy.MakeVertex(self.position[0],self.position[1], self.position[2] )
            circle1 = geompy.MakeCircle(circleCentre, None, self.radius)
            circle2 = geompy.MakeCircle(circleCentre, None, self.radius + self.filmThickness)
            face1 = geompy.MakeFaceWires([circle1, circle2], 1)
            vector1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
            circle3 = geompy.MakeCircle(circleCentre, vector1, 1)
            self.__ionomerPipe = geompy.MakePipe(face1, circle3)


from random import randint

particle = AgglomerateParticle( [randint(-5,+5),randint(-5,+5),randint(-5,+5),] , randint(1,5), 1, [randint(-5,+5),randint(-5,+5),randint(-5,+5)], randint(0,100))
particle.generateCarbonPlatinumGeometry()
particle.generateFilmGeometry()
geompy.addToStudy(particle.getCarbonParticle(), 'Particle')
geompy.addToStudy(particle.getIonomerPipe(), 'Film')

"""from time import time 

start = time()

filmOrient = [(math.pi)/2, 0, 0]
radius = 1
position = [0, 0 , 0]
filmThickness = 1

theta = float(math.pi/5)
zetha = theta/2 +filmOrient[0] -(math.pi)/2
gamma = theta - zetha
print theta ,zetha, gamma
arc1Pt1 = [-math.sin(zetha)*radius + position[0], math.cos(zetha)*radius + position[1], position[2]] 
arc1Pt2 = [math.sin(gamma)*radius  + position[0], math.cos(gamma)*radius  + position[1], position[2] ]

arc2Pt1 = [arc1Pt1[0] -math.sin(zetha)*filmThickness, arc1Pt1[1] +math.cos(zetha)*filmThickness, arc1Pt1[2] ]
arc2Pt2 = [arc1Pt2[0] +math.sin(gamma)*filmThickness, arc1Pt2[1] +math.cos(gamma)*filmThickness, arc1Pt2[2] ]

print arc1Pt1, arc1Pt2
print arc2Pt1, arc2Pt2
print time() - start"""