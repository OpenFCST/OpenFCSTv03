import smesh, geompy, SMESH
import numpy
import SALOMEDS
 
'''
Created on 3 Dec 2012

@author: Phil Wardlaw
'''


class mkSalomeMeshBase:
    '''
    classdocs
    '''

    

    def __init__(self, dim, gridDensity):
        '''
        Constructor
        '''
        
        self.subMeshes = []
        #self.compoundMesh = []
        self.__dimension = 2
        self.meshDensity = gridDensity #wire discretization length
        self.dimension = dim
        
    
    def generate(self):
        print "Generate function not yet implemented."
        raise Exception
    
    def getMesh(self):
        return self.compoundMesh
    
    def makeRectanglarMesh(self, dList,width, height):

        if len(dList) <3:
            dList.append(0)
            
        Vertex_1 = geompy.MakeVertex(dList[0], dList[1], dList[2])
        Vertex_2 = geompy.MakeVertex(dList[0] + width, dList[1], dList[2])
        Vertex_3 = geompy.MakeVertex(dList[0] + width, dList[1] + height, dList[2])
        Vertex_4 = geompy.MakeVertex(dList[0], dList[1] + height, dList[2])
        rect = geompy.MakeQuad4Vertices(Vertex_1, Vertex_2, Vertex_3, Vertex_4)
        Mesh_1 = smesh.Mesh(rect)
        Regular_1D = Mesh_1.Segment()
        Local_Length_1 = Regular_1D.LocalLength(self.meshDensity)
        Local_Length_1.SetPrecision( 1e-07 )
        Mesh_1.Quadrangle()
        Mesh_1.Compute()
        
        return Mesh_1 
    
    
    def mkCompound(self):
        if not self.subMeshes:
            print "Input mesh list is empty: Cannot build compound"
        else:
            self.compoundMesh=smesh.Concatenate([self.subMeshes[0].GetMesh()],1,1,1e-05)
            for i in range(len(self.subMeshes)-1):
            
                self.compoundMesh=smesh.Concatenate([self.compoundMesh.GetMesh() ,self.subMeshes[i+1].GetMesh()],1,1,1e-05) #check first one later
            #smesh.SetName(self.compoundMesh.GetMesh(),'Compound_Mesh')
            
    def delInternalEdges(self):
        
        if not self.compoundMesh:
            print "Cannot remove internal edges, mesh not initialized!"
            raise Exception
        
        
        search_filter = smesh.GetFilter(smesh.EDGE, smesh.FT_FreeBorders)
        external_borders = self.compoundMesh.GetIdsFromFilter(search_filter)
        
        #get all edges
        all_borders = self.compoundMesh.GetElementsByType(SMESH.EDGE)
        
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
        self.compoundMesh.RemoveElements(borders_to_remove)
        
        
    def delInternalFaces(self):
        
        if self.compoundMesh is None:
            print "Cannot remove internal faces, mesh not initialized!"
        
        else:
        
            search_filter = smesh.GetFilter(smesh.FACE, smesh.FT_FreeFaces)
            external_Faces = self.compoundMeshThreeD.GetIdsFromFilter(search_filter)
            
            #get all faces
            all_Faces = self.compoundMeshThreeD.GetElementsByType(SMESH.FACE)
            
            faces_to_remove = []
            
            #find the difference
            for f in all_Faces:
                if f in external_Faces:
                    pass
                else:
                    faces_to_remove.append(f)
            
            print "Removing internal faces:"
            print external_Faces
            print all_Faces
            print faces_to_remove
            
            #remove the faces
            self.compoundMeshThreeD.RemoveElements(faces_to_remove)
            
