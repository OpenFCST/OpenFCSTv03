'''
Python script for removing internal edge numbering in Salome
Created by: Philip Wardlaw 30 Oct 2012
'''


import geompy, salome, smesh, SMESH 


#Rename mesh_name to the desired mesh that you wish to remove internal edges from
#------------------------------------#
mesh_name = "Compound_Mesh_5"
#------------------------------------#

#Get mesh
mesh1= smesh.Mesh(salome.myStudy.FindObject(mesh_name).GetObject()) 

#Get external edges
search_filter = smesh.GetFilter(smesh.EDGE, smesh.FT_FreeBorders)
external_borders = mesh1.GetIdsFromFilter(search_filter)

#get all edges
all_borders = mesh1.GetElementsByType(SMESH.EDGE)

borders_to_remove = []

#find the difference
for b in all_borders:
  if b in external_borders:
    pass
  else:
    borders_to_remove.append(b)

print "Removing internal boarders:"
print borders_to_remove

#remove the edges
mesh1.RemoveElements(borders_to_remove)

#Refresh GUI
salome.sg.updateObjBrowser(1)


print "Done"
