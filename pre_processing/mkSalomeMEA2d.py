import smesh, geompy, SMESH
import numpy
import SALOMEDS
import sys
import os, inspect

print "="*50
print "= Script to create a 2D MEA with Salome"
print "-"*50

print "= - Set directory to find fcst Python Files"
fcst_dir = os.getenv('FCST_DIR')
if fcst_dir:
    print "=   Found FCST here: ", fcst_dir
    sys.path.append(fcst_dir + "/pre_processing")
else:
    print "=   FCST not found, inspect to guess path"
    print "=   do source PATH_TO_FCST/fcst_env.sh to set the approporate paths"
    print inspect.getfile(inspect.currentframe()) # script filename (usually with path)
    print os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    dirname= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    sys.path.append(dirname)

try:
    from mkSalomeMeshBase import mkSalomeMeshBase
except ImportError:
    print "FAIL!!!"
    sys.exit(1) 

#reload(mkSalomeMeshBase)

class mkSalomeMEA2d(mkSalomeMeshBase):
    
    def __init__(self,dim, gridDensity, anodebiPolarChannelWidth, anodebiPolarLandWidth , gdlWidth, mplWidth, aclWidth, membraneWidth, cclWidth, cathodebiPolarChannelWidth,cathodebiPolarLandWidth,totalHeight, anodechannelHeight,cathodechannelHeight):
        print "Child ", "-Constructor"
        mkSalomeMeshBase.__init__(self,dim, gridDensity)
        self.__widths = [anodebiPolarChannelWidth, anodebiPolarLandWidth, gdlWidth, mplWidth, aclWidth, membraneWidth, cclWidth, mplWidth, gdlWidth,cathodebiPolarLandWidth , cathodebiPolarChannelWidth ]
        self.__heights = [totalHeight, anodechannelHeight,cathodechannelHeight]
        if self.__widths[1]>0:       
            if self.__widths[0]>=self.__widths[1]:
                print "Incorrect input for anode channel width: assuming half of land width"
                self.__widths[0]==self.__widths[1]/2.0
        if self.__widths[-2]>0:
            if self.__widths[-1]>=self.__widths[-2]:
                print "Incorrect input for cathode channel width: assuming half of land width"
                self.__widths[-1]==self.__widths[-2]/2.0
        if self.__heights[0]<=0:
            print "Input Height for the MEA cannot be zero: Exiting code"
            raise Exception
        else:
            if self.__heights[0]<=self.__heights[1]: 
                print "Incorrect height for anode landing; Assuming half of total height as anode land height for anode"
                self.__heights[1]=self.__heights[0]/2
            if self.__heights[0]<=self.__heights[2]:
                print "Incorrect height for cathode landing; Assuming half of total height as cathode land height for anode"
                self.__heights[2]=self.__heights[0]/2
        print "Child", self.subMeshes
    
    def generate(self):
        if self.__widths[1]>0:
            for el in self.createAnodePlate():
                self.subMeshes.append(el)
        if self.__widths[-2]>0:
            for el in  self.createCathodePlate():
                self.subMeshes.append(el)
        
        
        
        x = self.__widths[1]
        materialIDCount = 3
        for w in self.__widths[2:-2]:
            
            rect = self.makeRectanglarMesh([x,0], w, self.__heights[0])
            
            
            a1 = rect.CreateEmptyGroup( SMESH.FACE, str(materialIDCount) )
            nbAdd = a1.AddFrom(rect.GetMesh() )
            a1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
            smesh.SetName(a1 ,  str(materialIDCount))
            self.subMeshes.append(rect)
            x += w
            materialIDCount +=1
            
        self.mkCompound()
        

    def createAnodePlate(self):
        #create the plate
        rect1 =  self.makeRectanglarMesh([0,0], self.__widths[1] -self.__widths[0], self.__heights[1])
        rect2 =  self.makeRectanglarMesh([0,self.__heights[1]], self.__widths[1], self.__heights[0]-self.__heights[1])
        
        #Need some more code:
        a1 = rect1.CreateEmptyGroup( SMESH.FACE, '1' )
        nbAdd = a1.AddFrom(rect1.GetMesh() )
        a1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
        a2 = rect2.CreateEmptyGroup( SMESH.FACE, '1' )
        nbAdd = a2.AddFrom(rect2.GetMesh() )
        a2.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
        
        
        smesh.SetName(a1, str(1))
        smesh.SetName(a2, str(1))
               
        #create the void
        rect3 =  self.makeRectanglarMesh([self.__widths[1] -self.__widths[0], 0], self.__widths[0], self.__heights[1])
        
        a3 = rect3.CreateEmptyGroup( SMESH.FACE, '2' )
        nbAdd = a3.AddFrom(rect3.GetMesh() )
        a3.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
        
        
        smesh.SetName(a3, str(2))
        
        return [rect1, rect2, rect3]
        
    def createCathodePlate(self):
        
        
        
        rightCornerX = sum(self.__widths) - self.__widths[0] - self.__widths[-1]
        
        
        #create the plate
        print "right corner",  rightCornerX
        print self.__widths[-2], self.__widths[-1]
        #print "op", [rightCornerX[0] ,0], -(self.__widths[1] -self.__widths[0]), self.__heights[0]
        rect1 =  self.makeRectanglarMesh([rightCornerX  - (self.__widths[-2] - self.__widths[-1]), 0 ], self.__widths[-2] - self.__widths[-1] , self.__heights[2])
        rect2 =  self.makeRectanglarMesh([rightCornerX - self.__widths[-2]  , self.__heights[2] ], self.__widths[-2], self.__heights[0]-self.__heights[2])
        
        a1 = rect1.CreateEmptyGroup( SMESH.FACE, '11' )
        nbAdd = a1.AddFrom(rect1.GetMesh() )
        a1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
        a2 = rect2.CreateEmptyGroup( SMESH.FACE, '11' )
        nbAdd = a2.AddFrom(rect2.GetMesh() )
        a2.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
        
        
        smesh.SetName(a1, str(1))
        smesh.SetName(a2, str(1))
               
        #create the void
        rect3 =  self.makeRectanglarMesh([rightCornerX -self.__widths[-2], 0], self.__widths[-1], self.__heights[2])
        a3 = rect3.CreateEmptyGroup( SMESH.FACE, '10' )
        nbAdd = a3.AddFrom(rect3.GetMesh() )
        a3.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
        
        
        smesh.SetName(a3, str(2))
        
        return [rect1, rect2, rect3]


if __name__ == '__main__':
    #dim, anodebiPolarChannelWidth, anodebiPolarLandWidth , gdlWidth, mplWidth, aclWidth, membraneWidth, cclWidth, cathodebiPolarChannelWidth,cathodebiPolarLandWidth,totalHeight, anodechannelHeight,cathodechannelHeight
    mea = mkSalomeMEA2d(2, 0.5, 1, 3, 3, 1, 1, 1, 1, 1, 2, 8, 3, 3)
    mea.generate()
    
    mea.delInternalEdges()
    mea.getMesh().ExportUNV('test.unv')

        
    
        