'''
Python script for removing internal edge numbering in Salome
Created by: Philip Wardlaw 30 Oct 2012
'''


import geompy, salome, smesh, SMESH 


#Rename mesh_name to the desired mesh that you wish to remove internal edges from
#------------------------------------#
mesh_name = "one_by_one_square_mesh.unv"
#------------------------------------#

#Get mesh
mesh1= smesh.Mesh(salome.myStudy.FindObject(mesh_name).GetObject()) 

#Get external edges

search_filter = smesh.GetFilter(smesh.FACE, smesh.FT_FreeFaces)
external_faces = mesh1.GetIdsFromFilter(search_filter)

#get all faces
all_Faces = mesh1.GetElementsByType(SMESH.FACE)

faces_to_remove = []

#find the difference
for f in all_Faces:
    if f in external_faces:
        pass
    else:
        faces_to_remove.append(f)

print "Removing internal faces:"
print faces_to_remove


#remove the edges
mesh1.RemoveElements(faces_to_remove)

#Refresh GUI
salome.sg.updateObjBrowser(1)


print "Done"
