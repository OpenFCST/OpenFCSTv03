'''
Created on 4 Dec 2012

@author: Phil Wardlaw
'''

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
    from mkSalomeMEA2d import mkSalomeMEA2d
except ImportError:
    print "FAIL!!!"
    sys.exit(1) 

class mkSalomeMEA3d(mkSalomeMEA2d):
    '''
    classdocs
    '''


    def __init__(self,dim, gridDensity, anodebiPolarChannelWidth, anodebiPolarLandWidth, gdlWidth, mplWidth, aclWidth, membraneWidth, cclWidth, cathodebiPolarChannelWidth, cathodebiPolarLandWidth, totalHeight, anodechannelHeight, cathodechannelHeight,  depth):
        mkSalomeMEA2d.__init__(self, dim, gridDensity, anodebiPolarChannelWidth, anodebiPolarLandWidth, gdlWidth, mplWidth, aclWidth, membraneWidth, cclWidth, cathodebiPolarChannelWidth, cathodebiPolarLandWidth, totalHeight, anodechannelHeight, cathodechannelHeight)
        self.__depth = depth
        self.compoundMeshThreeD = smesh.Mesh()
        
    def generate(self):
        mkSalomeMEA2d.generate(self)
        if self.__depth % self.meshDensity:
            print "Warning: z depth is not equally divided by mesh density (Try make the depth an multiple of the mesh density)"
        self.addToCompoundThreeDMesh(self.compoundMesh.ExtrusionSweepObject2D( self.compoundMesh, SMESH.DirStruct( SMESH.PointStruct ( 0, 0, self.meshDensity )), int(self.__depth/self.meshDensity) ,True))
    
    
    def addToCompoundThreeDMesh(self, subThreeDMeshes):
            for mesh in subThreeDMeshes:
                self.compoundMeshThreeD = smesh.Concatenate([self.compoundMeshThreeD.GetMesh() ,mesh.GetMesh()],1,1,1e-05)
                
    def getMesh(self):
        return self.compoundMeshThreeD
    
    #overide base function
    def delInternalEdges(self):
        search_filter = smesh.GetFilter(smesh.EDGE, smesh.FT_FreeBorders)
        external_borders = self.compoundMeshThreeD.GetIdsFromFilter(search_filter)
        
        #get all edges
        all_borders = self.compoundMeshThreeD.GetElementsByType(SMESH.EDGE)
        
        borders_to_remove = []
        
        #find the difference
        for b in all_borders:
            if b in external_borders:
                pass
            else:
                borders_to_remove.append(b)
        
        print "Removing internal edges:"
        print borders_to_remove
        
        #remove the edges
        self.compoundMeshThreeD.RemoveElements(borders_to_remove)

if __name__ == '__main__':
    #dim, anodebiPolarChannelWidth, anodebiPolarLandWidth , gdlWidth, mplWidth, aclWidth, membraneWidth, cclWidth, cathodebiPolarChannelWidth,cathodebiPolarLandWidth,totalHeight, anodechannelHeight,cathodechannelHeight
    mea = mkSalomeMEA3d(2, 0.125, 1, 3, 3, 1, 1, 1, 1, 1, 2, 8, 3, 3, 0.375)
    mea.generate()
    
    mea.delInternalEdges()
    mea.delInternalFaces()
    mea.getMesh().ExportUNV('test.unv')    