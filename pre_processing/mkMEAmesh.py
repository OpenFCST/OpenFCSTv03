

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy



####################################################
##       Begin of variables section      ##
####################################################
width = 0.06
C_plate_h = 0.004
A_plate_h = 0.003
C_GDL_h = 0.001
A_GDL_h = 0.001
C_MPL_h = 0.002
A_MPL_h = 0.002
C_CL_h = 0.001
A_CL_h = 0.001
Mem_h = 0.002
C_Ch_h = 0.001
C_Ch_w = 0.002
A_Ch_h = 0.001
A_Ch_w = 0.002
A_Cool_h = 0.001
A_Cool_w = 0.002
C_land_h = 0.001
C_land_w = width-C_Ch_w
A_land_h = 0.0001
A_land_w = width-A_Ch_w
mesh_length=0.001
####################################################
##        End of variables section       ##
####################################################


###
### GEOM component
###

import GEOM
import geompy
import math
import SALOMEDS


geompy.init_geom(theStudy)

Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(0, width, 0)
Vertex_3 = geompy.MakeVertex(C_GDL_h, width, 0)
Vertex_4 = geompy.MakeVertex(C_GDL_h, 0, 0)
Vertex_5 = geompy.MakeVertex(C_GDL_h+C_MPL_h, 0, 0)
Vertex_6 = geompy.MakeVertex(C_GDL_h+C_MPL_h, width, 0)
Vertex_7 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h, width, 0)
Vertex_8 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h, 0, 0)
Vertex_9 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h, 0, 0)
Vertex_10 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h, width, 0)
Vertex_11 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h+A_CL_h, width, 0)
Vertex_12 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h+A_CL_h, 0, 0)
Vertex_13 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h+A_CL_h+A_MPL_h, 0, 0)
Vertex_14 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h+A_CL_h+A_MPL_h, width, 0)
Vertex_15 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h+A_CL_h+A_MPL_h+A_GDL_h, width, 0)
Vertex_16 = geompy.MakeVertex(C_GDL_h+C_MPL_h+C_CL_h+Mem_h+A_CL_h+A_MPL_h+A_GDL_h, 0, 0)
Quadrangle_Face_1 = geompy.MakeQuad4Vertices(Vertex_1, Vertex_2, Vertex_3, Vertex_4)
Quadrangle_Face_2 = geompy.MakeQuad4Vertices(Vertex_3, Vertex_4, Vertex_5, Vertex_6)
Quadrangle_Face_3 = geompy.MakeQuad4Vertices(Vertex_5, Vertex_6, Vertex_7, Vertex_8)
Quadrangle_Face_4 = geompy.MakeQuad4Vertices(Vertex_7, Vertex_8, Vertex_9, Vertex_10)
Quadrangle_Face_5 = geompy.MakeQuad4Vertices(Vertex_9, Vertex_10, Vertex_11, Vertex_12)
Quadrangle_Face_6 = geompy.MakeQuad4Vertices(Vertex_11, Vertex_12, Vertex_13, Vertex_14)
Quadrangle_Face_7 = geompy.MakeQuad4Vertices(Vertex_13, Vertex_14, Vertex_15, Vertex_16)
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Vertex_10, 'Vertex_10' )
geompy.addToStudy( Vertex_11, 'Vertex_11' )
geompy.addToStudy( Vertex_12, 'Vertex_12' )
geompy.addToStudy( Vertex_13, 'Vertex_13' )
geompy.addToStudy( Vertex_14, 'Vertex_14' )
geompy.addToStudy( Vertex_15, 'Vertex_15' )
geompy.addToStudy( Vertex_16, 'Vertex_16' )
geompy.addToStudy( Quadrangle_Face_1, 'Quadrangle Face_1' )
geompy.addToStudy( Quadrangle_Face_2, 'Quadrangle Face_2' )
geompy.addToStudy( Quadrangle_Face_3, 'Quadrangle Face_3' )
geompy.addToStudy( Quadrangle_Face_4, 'Quadrangle Face_4' )
geompy.addToStudy( Quadrangle_Face_5, 'Quadrangle Face_5' )
geompy.addToStudy( Quadrangle_Face_6, 'Quadrangle Face_6' )
geompy.addToStudy( Quadrangle_Face_7, 'Quadrangle Face_7' )



###
### SMESH component
###

import smesh, SMESH, SALOMEDS

aMeasurements = smesh.CreateMeasurements()
smesh.SetCurrentStudy(theStudy)
import StdMeshers
Mesh_1 = smesh.Mesh(Quadrangle_Face_1)
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(mesh_length)
Local_Length_1.SetPrecision( 1e-07 )
Quadrangle_2D = Mesh_1.Quadrangle()
isDone = Mesh_1.Compute()
Mesh_2 = smesh.Mesh(Quadrangle_Face_2)
status = Mesh_2.AddHypothesis(Local_Length_1)
Regular_1D_1 = Mesh_2.Segment()
Quadrangle_2D_1 = Mesh_2.Quadrangle()
isDone = Mesh_2.Compute()
Mesh_3 = smesh.Mesh(Quadrangle_Face_3)
status = Mesh_3.AddHypothesis(Local_Length_1)
Regular_1D_2 = Mesh_3.Segment()
Quadrangle_2D_2 = Mesh_3.Quadrangle()
isDone = Mesh_3.Compute()
Mesh_4 = smesh.Mesh(Quadrangle_Face_4)
status = Mesh_4.AddHypothesis(Local_Length_1)
Regular_1D_3 = Mesh_4.Segment()
Quadrangle_2D_3 = Mesh_4.Quadrangle()
isDone = Mesh_4.Compute()
Mesh_5 = smesh.Mesh(Quadrangle_Face_5)
status = Mesh_5.AddHypothesis(Local_Length_1)
Regular_1D_4 = Mesh_5.Segment()
Quadrangle_2D_4 = Mesh_5.Quadrangle()
isDone = Mesh_5.Compute()
Mesh_6 = smesh.Mesh(Quadrangle_Face_6)
status = Mesh_6.AddHypothesis(Local_Length_1)
Regular_1D_5 = Mesh_6.Segment()
Quadrangle_2D_5 = Mesh_6.Quadrangle()
isDone = Mesh_6.Compute()
Mesh_7 = smesh.Mesh(Quadrangle_Face_7)
status = Mesh_7.AddHypothesis(Local_Length_1)
Regular_1D_6 = Mesh_7.Segment()
Quadrangle_2D_6 = Mesh_7.Quadrangle()
isDone = Mesh_7.Compute()
a1 = Mesh_1.CreateEmptyGroup( SMESH.FACE, '1' )
nbAdd = a1.AddFrom( Mesh_1.GetMesh() )
a1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a2 = Mesh_2.CreateEmptyGroup( SMESH.FACE, '2' )
nbAdd = a2.AddFrom( Mesh_2.GetMesh() )
a2.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a3 = Mesh_3.CreateEmptyGroup( SMESH.FACE, '3' )
nbAdd = a3.AddFrom( Mesh_3.GetMesh() )
a3.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a4 = Mesh_4.CreateEmptyGroup( SMESH.FACE, '4' )
nbAdd = a4.AddFrom( Mesh_4.GetMesh() )
a4.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a5 = Mesh_5.CreateEmptyGroup( SMESH.FACE, '5' )
nbAdd = a5.AddFrom( Mesh_5.GetMesh() )
a5.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a6 = Mesh_6.CreateEmptyGroup( SMESH.FACE, '6' )
nbAdd = a6.AddFrom( Mesh_6.GetMesh() )
a6.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a7 = Mesh_7.CreateEmptyGroup( SMESH.FACE, '7' )
nbAdd = a7.AddFrom( Mesh_7.GetMesh() )
a7.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a1_1 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, '1' )
nbAdd = a1_1.Add( [ 1,2,3,4,5,6,7,8,9,10,11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30] )
a1_1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a2_1 = Mesh_1.CreateEmptyGroup( SMESH.EDGE, '2' )
nbAdd = a2_1.Add( [ 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60 ] )
a2_1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a3_1 = Mesh_7.CreateEmptyGroup( SMESH.EDGE, '3' )
nbAdd = a3_1.Add( [65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95] )
a3_1.SetColor( SALOMEDS.Color( 0, 0.666667, 1 ))
a4_1 = Mesh_7.CreateEmptyGroup( SMESH.EDGE, '4' )
nbAdd = a4_1.Add( [ 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124] )

[ a1, a1_1, a2_1 ] = Mesh_1.GetGroups()
[ a2 ] = Mesh_2.GetGroups()
[ a3 ] = Mesh_3.GetGroups()
[ a4 ] = Mesh_4.GetGroups()
[ a5 ] = Mesh_5.GetGroups()
[ a6 ] = Mesh_6.GetGroups()
[ a7, a3_1, a4_1 ] = Mesh_7.GetGroups()
Compound_Mesh_1 = smesh.Concatenate([Mesh_1.GetMesh(), Mesh_2.GetMesh(), Mesh_3.GetMesh(), Mesh_4.GetMesh(), Mesh_5.GetMesh(), Mesh_6.GetMesh(), Mesh_7.GetMesh()], 1, 1, 1e-05)
[ a1_2, a1_3, a2_2, a2_3, a3_2, a4_2, a5_1, a6_1, a7_1, a3_3, a4_3 ] = Compound_Mesh_1.GetGroups()

# isDone = Compound_Mesh_1.RemoveElements( [ 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124 ] )

## set object names
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Mesh_2.GetMesh(), 'Mesh_2')
smesh.SetName(Mesh_3.GetMesh(), 'Mesh_3')
smesh.SetName(Mesh_4.GetMesh(), 'Mesh_4')
smesh.SetName(Mesh_5.GetMesh(), 'Mesh_5')
smesh.SetName(Mesh_6.GetMesh(), 'Mesh_6')
smesh.SetName(Mesh_7.GetMesh(), 'Mesh_7')
smesh.SetName(a1, '1')
smesh.SetName(a2, '2')
smesh.SetName(a3, '3')
smesh.SetName(a4, '4')
smesh.SetName(a5, '5')
smesh.SetName(a6, '6')
smesh.SetName(a7, '7')
smesh.SetName(a1_1, '1')
smesh.SetName(a2_1, '2')
smesh.SetName(a3_1, '3')
smesh.SetName(a4_1, '4')
smesh.SetName(Compound_Mesh_1.GetMesh(), 'Compound_Mesh_1')
smesh.SetName(a1_2, '1')
smesh.SetName(a1_3, '1')
smesh.SetName(a2_2, '2')
smesh.SetName(a2_3, '2')
smesh.SetName(a3_2, '3')
smesh.SetName(a4_2, '4')
smesh.SetName(a5_1, '5')
smesh.SetName(a6_1, '6')
smesh.SetName(a7_1, '7')
smesh.SetName(a3_3, '3')
smesh.SetName(a4_3, '4')

## convert mesh to unv format##

Compound_Mesh_1.RenumberNodes()
Compound_Mesh_1.RenumberElements()
Compound_Mesh_1.ExportUNV( r'/global/home/tinsho/fcst/data/mea/parametric/salome_mesh_final.unv' )
