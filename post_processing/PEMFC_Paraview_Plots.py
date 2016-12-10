"""
Simple Paraview script for creating batchs of contour plots.
Work in progress, subject to a lot of change over November 2013

Works for vtk output from PEMFC app.

To use edit the following 4 variables:

fileName: the vtk to be proccessed
exportDir: the location where the images will be placed (make sure the folder already exists)
plotName: just a string used to name the files
layer: 4.0 is the catalyst (the only layer this was developed for)

Your milage may vary. 

Tips: 
    Before you run this script go edit>delete all **(Important otherwise ranges will be off!)**
    To run the script open Paraview>tools>python shell>run script>select this file

"""

#---------------------------------------------------------------------------------------------------------
# If using all the user needs to do is follow the steps above.
#=========================================================================================================
exportDir = '/home/wardlawp/Code_workspace/FCST/data/AgglomeratePaper/sigma_P/3_two_orders_smaller_by5/post/'
fileName = exportDir +'1.vtk'
plotName = "Case1_1.0A"



# Is the same as above however during run time you pass the name of the vtk file (Voltage 0.000000, cycle number, Name of file you want to save plots under.).
#=========================================================================================================
#voltage = raw_input("Voltage plz:")
#cycleN = raw_input("Cycle number plz:")
#Name = raw_input("Current [A] plz:")
#exportDir = '/home/kdomican/FCST_2013_09SEP_18/trunk/data/mea/Kailyn/MEA_Param_Homo_MPL_DTK_CD_125_Fitted_th_RHsecond/MEA_40_Nafion/Anode_DP_Agg__Cathode_DT_Agg/'
#fileName = exportDir +'fuel_cell_solution_DataFile_000'+cycleN+'_Cycle_1__V_cell_-'+voltage+'.vtk'
#plotName = ""+Name+"A"

#---------------------------------------------------------------------------------------------------------


layer = 4.0


def cleanEverything():
    for key in GetSources():
        s = FindSource(key[0])
        Delete(s)
    
    for key in GetRepresentations():
        r = FindSource(key[0])
        Delete(r)



try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

cleanEverything()




fuel_cell_solution_DataFile_00017_Cycle_1__V_cell_0_548000_vtk = LegacyVTKReader( FileNames=[fileName] )


RenderView1 = GetRenderView()
RenderView1.InteractionMode = '2D'

RenderView1.Background = [1.0, 1.0, 1.0]
RenderView1.CenterOfRotation = [0.8333250000000001, 0.05, 0.0]
RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
RenderView1.OrientationAxesVisibility = 0


#---------------------------------------------------------------------------------------------------------
#Conventional CL
#==============================================================
RenderView1.CameraFocalPoint = [0.8333250000000001, 0.05, 0.0]
RenderView1.CameraParallelScale = 0.051538820320220766
RenderView1.CameraPosition = [0.8333250000000001, 0.05, 0.19913071041509228]
RenderView1.CameraClippingRange = [0.19713940331094135, 0.20211767107131867]

# Low Loading CL
#==============================================================
#RenderView1.CameraFocalPoint = [4.91325, 0.05, 0.0]
#RenderView1.CameraParallelScale = 0.050717354031928785
#RenderView1.CameraPosition = [4.91325, 0.05, 0.19936718414350557]
#RenderView1.CameraClippingRange = [0.1973735123020705, 0.20235769190565817]

#---------------------------------------------------------------------------------------------------------


Calculator5 = Calculator()
Calculator5.Function = '-1.16151-protonic_electrical_potential+electronic_electrical_potential'
Calculator5.ResultArrayName = 'Overpotential'


Threshold1 = Threshold()
Threshold1.Scalars = ['POINTS', 'cell_material_ids']
Threshold1.ThresholdRange = [layer, layer]

DataRepresentation2 = Show()
DataRepresentation2.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation2.AmbientColor = [0.0, 0.0, 0.0]



#---------------------------------------------------------------------------------------------------------
#Conventional CL
#==============================================================
#DataRepresentation2.Scale = [150.0, 1.0, 1.0]


# Low Loading CL
#==============================================================
DataRepresentation2.Scale = [25.0, 1.0, 1.0]

#---------------------------------------------------------------------------------------------------------





#######################################protonic_electrical_potential

rangePhiM = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation("protonic_electrical_potential").GetComponentRange(0)
a1_protonic_electrical_potential_PVLookupTable = GetLookupTableForArray( "protonic_electrical_potential", 1, RGBPoints=[rangePhiM[0], 0.0, 0.0, 1.0, rangePhiM[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )
a1_protonic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangePhiM[0], 0.0, 0.5, 0.0, rangePhiM[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation2 = CreateScalarBar( Title='Phi_M', Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_protonic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation2)



DataRepresentation2.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'protonic_electrical_potential'
DataRepresentation2.LookupTable = a1_protonic_electrical_potential_PVLookupTable

a1_protonic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName+ '_phiM.png')

ScalarBarWidgetRepresentation2.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation2.Enabled = 0
ScalarBarWidgetRepresentation2.Visibility = 0
ScalarBarWidgetRepresentation2.LabelColor = [0.0, 0.0, 0.0]

#######################################electronic_electrical_potential
rangePhiS = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation("electronic_electrical_potential").GetComponentRange(0)
a1_electronic_electrical_potential_PVLookupTable = GetLookupTableForArray( "electronic_electrical_potential", 1, RGBPoints=[rangePhiS[0], 0.0, 0.0, 1.0, rangePhiS[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_electronic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangePhiS[0], 0.0, 0.5, 0.0, rangePhiS[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation3 = CreateScalarBar( Title='Phi_S', Position2=[0.13, 0.49999999999999956], Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_electronic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0],  Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation3)


DataRepresentation2.ScalarOpacityFunction = a1_electronic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'electronic_electrical_potential'

DataRepresentation2.LookupTable = a1_electronic_electrical_potential_PVLookupTable

a1_electronic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_electronic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName +'_phiS.png')

ScalarBarWidgetRepresentation3.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation3.Enabled = 0
ScalarBarWidgetRepresentation3.Visibility = 0
ScalarBarWidgetRepresentation3.LabelColor = [0.0, 0.0, 0.0]

#######################################Agglomerate Effectiveness

rangeEff = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation("agglomerate_effectiveness").GetComponentRange(0)
a1_agglomerate_effectiveness_PVLookupTable = GetLookupTableForArray( "agglomerate_effectiveness", 1, RGBPoints=[rangeEff[0], 0.0, 0.0, 1.0, rangeEff[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_agglomerate_effectiveness_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangeEff[0], 0.0, 0.5, 0.0, rangeEff[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation5 = CreateScalarBar( Title='Agg Eff [%]', Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_agglomerate_effectiveness_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation5)

DataRepresentation2.ScalarOpacityFunction = a1_agglomerate_effectiveness_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'agglomerate_effectiveness'
DataRepresentation2.LookupTable = a1_agglomerate_effectiveness_PVLookupTable


a1_agglomerate_effectiveness_PVLookupTable.ScalarOpacityFunction = a1_agglomerate_effectiveness_PiecewiseFunction

WriteImage(exportDir + plotName + '_eff.png')


ScalarBarWidgetRepresentation5.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation5.Enabled = 0
ScalarBarWidgetRepresentation5.Visibility = 0
ScalarBarWidgetRepresentation5.LabelColor = [0.0, 0.0, 0.0]

#######################################Membrane Water Content

scale1 = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation( "membrane_water_content").GetComponentRange(0)
a1_membrane_water_content_PVLookupTable = GetLookupTableForArray( "membrane_water_content", 1, RGBPoints=[scale1[0], 0.0, 0.0, 1.0, scale1[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )
a1_membrane_water_content_PiecewiseFunction = CreatePiecewiseFunction( Points=[scale1[0], 0.0, 0.5, 0.0, scale1[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation4 = CreateScalarBar( Title='   Water Content', Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_membrane_water_content_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation4)


DataRepresentation2.ScalarOpacityFunction = a1_membrane_water_content_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'membrane_water_content'
DataRepresentation2.LookupTable = a1_membrane_water_content_PVLookupTable


a1_membrane_water_content_PVLookupTable.ScalarOpacityFunction = a1_membrane_water_content_PiecewiseFunction

WriteImage(exportDir + plotName+ '_water.png')

ScalarBarWidgetRepresentation4.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation4.Enabled = 0
ScalarBarWidgetRepresentation4.Visibility = 0
ScalarBarWidgetRepresentation4.LabelColor = [0.0, 0.0, 0.0]

#######################################Relative Humidity

rangeRelH = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation( "relative_humidity").GetComponentRange(0)
a1_protonic_electrical_potential_PVLookupTable = GetLookupTableForArray( "relative_humidity", 1, RGBPoints=[rangeRelH[0], 0.0, 0.0, 1.0, rangeRelH[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )
a1_protonic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangeRelH[0], 0.0, 0.5, 0.0, rangeRelH[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation6 = CreateScalarBar( Title='      Relative Humidity', Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_protonic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation6)


DataRepresentation2.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'relative_humidity'
DataRepresentation2.LookupTable = a1_protonic_electrical_potential_PVLookupTable

a1_protonic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName+ '_relHumidity.png')


ScalarBarWidgetRepresentation6.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation6.Enabled = 0
ScalarBarWidgetRepresentation6.Visibility = 0
ScalarBarWidgetRepresentation6.LabelColor = [0.0, 0.0, 0.0]

#######################################water_molar_fraction
rangeWmF = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation( "water_molar_fraction").GetComponentRange(0)
a1_protonic_electrical_potential_PVLookupTable = GetLookupTableForArray( "water_molar_fraction", 1, RGBPoints=[rangeWmF[0], 0.0, 0.0, 1.0, rangeWmF[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_protonic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangeWmF[0], 0.0, 0.5, 0.0, rangeWmF[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation7 = CreateScalarBar( Title='         H2O Molar Fraction', Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_protonic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation7)


DataRepresentation2.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'water_molar_fraction'
DataRepresentation2.LookupTable = a1_protonic_electrical_potential_PVLookupTable


a1_protonic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName+ '_waterMolarFraction.png')


ScalarBarWidgetRepresentation7.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation7.Enabled = 0
ScalarBarWidgetRepresentation7.Visibility = 0
ScalarBarWidgetRepresentation7.LabelColor = [0.0, 0.0, 0.0]
#######################################"oxygen_molar_fraction"
rangeO2 = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation( "oxygen_molar_fraction").GetComponentRange(0)
a1_protonic_electrical_potential_PVLookupTable = GetLookupTableForArray("oxygen_molar_fraction", 1, RGBPoints=[rangeO2[0], 0.0, 0.0, 1.0, rangeO2[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_protonic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangeO2[0], 0.0, 0.5, 0.0, rangeO2[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation8 = CreateScalarBar( Title="         O2 Molar Fraction", Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_protonic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation8)


DataRepresentation2.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = "oxygen_molar_fraction"
DataRepresentation2.LookupTable = a1_protonic_electrical_potential_PVLookupTable

a1_protonic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName+ '_o2.png')

ScalarBarWidgetRepresentation8.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation8.Enabled = 0
ScalarBarWidgetRepresentation8.Visibility = 0
ScalarBarWidgetRepresentation8.LabelColor = [0.0, 0.0, 0.0]

####################################### Overpotential & Contour lines & Auto Cropping





rangOverPot = Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation('Overpotential').GetComponentRange(0)
a1_protonic_electrical_potential_PVLookupTable = GetLookupTableForArray("Overpotential", 1, RGBPoints=[rangOverPot[0], 0.0, 0.0, 1.0, rangOverPot[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_protonic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangOverPot[0], 0.0, 0.5, 0.0, rangOverPot[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation9 = CreateScalarBar( Title="Overpotential", Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_protonic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation9)


DataRepresentation2.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = "Overpotential"
DataRepresentation2.LookupTable = a1_protonic_electrical_potential_PVLookupTable

a1_protonic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName+ '_overPot.png')

ScalarBarWidgetRepresentation9.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation9.Enabled = 0
ScalarBarWidgetRepresentation9.Visibility = 0
ScalarBarWidgetRepresentation9.LabelColor = [0.0, 0.0, 0.0]


#################################

#######################################"current_density"
rangeI= Threshold1.GetDataInformation().GetPointDataInformation().GetArrayInformation( "current_density").GetComponentRange(0)
a1_protonic_electrical_potential_PVLookupTable = GetLookupTableForArray("current_density", 1, RGBPoints=[rangeI[0], 0.0, 0.0, 1.0, rangeI[1], 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.4980392156862745, 0.4980392156862745, 0.4980392156862745], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_protonic_electrical_potential_PiecewiseFunction = CreatePiecewiseFunction( Points=[rangeI[0], 0.0, 0.5, 0.0, rangeI[1], 1.0, 0.5, 0.0] )

ScalarBarWidgetRepresentation10 = CreateScalarBar( Title="   Current Density", Enabled=1, LabelFontSize=12, LabelColor=[0.0, 0.0, 0.0], LookupTable=a1_protonic_electrical_potential_PVLookupTable, TitleFontSize=12, TitleColor=[0.0, 0.0, 0.0], Position=[0.6, 0.21059113300492627] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation10)


DataRepresentation2.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction
DataRepresentation2.ColorArrayName = "current_density"
DataRepresentation2.LookupTable = a1_protonic_electrical_potential_PVLookupTable

a1_protonic_electrical_potential_PVLookupTable.ScalarOpacityFunction = a1_protonic_electrical_potential_PiecewiseFunction

WriteImage(exportDir + plotName+ '_i.png')

ScalarBarWidgetRepresentation10.TitleColor = [0.0, 0.0, 0.0]
ScalarBarWidgetRepresentation10.Enabled = 0
ScalarBarWidgetRepresentation10.Visibility = 0
ScalarBarWidgetRepresentation10.LabelColor = [0.0, 0.0, 0.0]
##################################################

Render()
