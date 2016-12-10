# -*- coding: utf-8 -*-

r"""
*************************************************************************
:mod:`PythonFCST.mesh.GridGenerator`: GridGenerator for VTK Meshes
*************************************************************************

.. module:: PythonFCST.mesh.GridGenerator

"""

import numpy as np
try:
    from tvtk.api import tvtk
except:
    try:
        from enthought.tvtk.api import tvtk
    except:
        print "TVTK package not found"


class PhaseGenerator:
    r"""
    This class generates a VTK mesh given a 3D numpy array.
    
    Parameters
    ----------
    image : ndarray
        3D image to be written to a VTK file
    
    filename : string
        Output filename for the VTK structure
       
    
    
    Examples
    --------
    
    >>> import PythonFCST as fcst
    >>> import fcst.mesh as mesh
    >>> crunch = mesh.GridGenerator(image,filename)
    >>> crunch.write()
    
    .. warning:: This documentation is not yet finished.
    
    .. todo:: 1. Finish documentation
        
    """
    def __init__(self,image, filename, scale=[1E-6,1e-6,1E-6],**kwargs):
        
        self.image=image
        xdim,ydim,zdim=np.shape(image)
        self.filename=filename
        x,y,z=np.mgrid[0:xdim+1,0:ydim+1,0:zdim+1]        
        self.x=x.flatten()
        self.y=y.flatten()
        self.z=z.flatten()
        self.scale=scale

        self._indent='\t'
        
    
    def generatePoints(self):

        points=np.array(np.dstack([self.x,self.y,self.z]),'f')
        points=points.reshape(self.x.shape[0],3)
        
        return points

    
    def generate3dCells(self,points):
        
        boundaryz =  np.max(self.z)
        boundaryx =  np.max(self.x)
        boundaryy =  np.max(self.y)
        x=self.x.reshape(self.x.shape[0],1)
        y=self.y.reshape(self.y.shape[0],1)
        z=self.z.reshape(self.z.shape[0],1) 
        
        
        a=np.max(self.z)+1
        b=(1+np.max(self.z))*(np.max(self.y)+1)

        count = self.x.shape[0]
        k = np.arange(count)
        i2=np.asarray([i for i in k if boundaryz not in z[i] and boundaryy not in y[i] and boundaryx not in x[i]])
        i1=np.ones(np.size(i2))*8
        i3=np.asarray([i+b for i in i2])
        i4=np.asarray([i+a+b for i in i2])
        i5=np.asarray([i+a for i in i2])
        i6=np.asarray([i+1 for i in i2])
        i7=np.asarray([i+b+1 for i in i2])
        i8=np.asarray([i+a+b+1 for i in i2])
        i9=np.asarray([i+a+1 for i in i2])
        
                
        cell_data=self.image
        cell_data=cell_data.flatten()
        
        index=np.where(cell_data>0)
        
        
        i1=np.delete(i1,index)
        i2=np.delete(i2,index)
        i3=np.delete(i3,index)
        i4=np.delete(i4,index)
        i5=np.delete(i5,index)
        i6=np.delete(i6,index)
        i7=np.delete(i7,index)
        i8=np.delete(i8,index)
        i9=np.delete(i9,index)
        cell_data=np.delete(cell_data,index)
        
        
        cellcount=i1.shape[0]
        
        cells=np.array(np.dstack([i1,i2,i3,i4,i5,i6,i7,i8,i9]))
        cells=cells.reshape(i1.shape[0],9)
        
        
        cells=cells.flatten()
        hex_type = tvtk.Hexahedron().cell_type # VTK_HEXAHEDRON == 12
        
        cell_types=np.ones(cellcount)*hex_type
        
        offset=np.arange(0,cells.shape[0],9)
        
        
        
        
        return [cell_types, offset, cellcount, cells, cell_data]
        
    def delete2dindex(self,i1,i2,i3,i4,i5,data,index):
        
        i1=np.delete(i1,index)
        i2=np.delete(i2,index)
        i3=np.delete(i3,index)
        i4=np.delete(i4,index)
        i5=np.delete(i5,index)
        data = np.delete(data,index)
        return i1,i2,i3,i4,i5,data
        
        
        
    def generate2dCells(self,points):
        
        boundaryz =  np.max(self.z)
        boundaryy =  np.max(self.y)
        boundaryx =  np.max(self.x)
        x=self.x.reshape(self.x.shape[0],1)
        y=self.y.reshape(self.y.shape[0],1)
        z=self.z.reshape(self.z.shape[0],1)        
        
        
        count = self.x.shape[0]
        k = np.arange(count)
        a=np.max(z)+1
        b=(1+np.max(z))*(np.max(y)+1)

        #Y-Z Plane boundaries
        
        
        ib2=np.asarray([i for i in k if 0 in x[i] and boundaryz not in z[i] and boundaryy not in y[i]])
        ib1=np.ones(np.size(ib2))*4
        ib3=np.asarray([i+1 for i in ib2])
        ib4=np.asarray([i+a+1 for i in ib2])
        ib5=np.asarray([i+a for i in ib2])
        data1 = self.image[0,:,:].copy()
        data1[:,:] = 7        
        data1[self.image[0,:,:]==0] = 3
        data1 = data1.flatten()
        index=np.where(data1==7)
        ib1,ib2,ib3,ib4,ib5,data1=self.delete2dindex(ib1,ib2,ib3,ib4,ib5,data1,index)
        
        #data1=np.ones(ib1.shape[0])
        
        ob2=np.asarray([i for i in k if boundaryx in x[i] and boundaryy not in y[i] and boundaryz not in z[i]])
        ob1=np.ones(np.size(ob2))*4
        ob3=np.asarray([i+1 for i in ob2])
        ob4=np.asarray([i+a+1 for i in ob2])
        ob5=np.asarray([i+a for i in ob2])
        data2 = self.image[-1,:,:].copy()
        data2[:,:] = 8        
        data2[self.image[-1,:,:]==0] = 4
        data2 = data2.flatten()
        index=np.where(data2==8)
        ob1,ob2,ob3,ob4,ob5,data2=self.delete2dindex(ob1,ob2,ob3,ob4,ob5,data2,index)
        #data2=np.ones(ob1.shape[0])*2
            
        
        innercellsX=np.array(np.dstack([ib1,ib2,ib3,ib4,ib5]))
        innercellsX=innercellsX.reshape(ib1.shape[0],5)
        
        outercellsX=np.array(np.dstack([ob1,ob2,ob3,ob4,ob5]))
        outercellsX=outercellsX.reshape(ob1.shape[0],5)
        
        #X-Z Plane boundaries
        
        count = self.y.shape[0]
        k = np.arange(count)
        
        
        ib2=np.asarray([i for i in k if 0 in y[i] and boundaryz not in z[i] and boundaryx not in x[i]])
        ib1=np.ones(np.size(ib2))*4
        ib3=np.asarray([i+1 for i in ib2])
        ib4=np.asarray([i+b+1 for i in ib2])
        ib5=np.asarray([i+b for i in ib2])
        data3 = self.image[:,0,:].copy()
        data3[:,:] = 9        
        data3[self.image[:,0,:]==0] = 1
        data3 = data3.flatten()
        index=np.where(data3==9)
        ib1,ib2,ib3,ib4,ib5,data3=self.delete2dindex(ib1,ib2,ib3,ib4,ib5,data3,index)
        #data3=np.ones(ob1.shape[0])*3
        
        ob2=np.asarray([i for i in k if boundaryy in y[i] and boundaryz not in z[i] and boundaryx not in x[i]])
        ob1=np.ones(np.size(ob2))*4
        ob3=np.asarray([i+1 for i in ob2])
        ob4=np.asarray([i+b+1 for i in ob2])
        ob5=np.asarray([i+b for i in ob2])
        data4 = self.image[:,-1,:].copy()
        data4[:,:] = 10        
        data4[self.image[:,-1,:]==0] = 2
        data4 = data4.flatten()
        #data4=np.ones(ob1.shape[0])*4      
        index=np.where(data4==10)
        ob1,ob2,ob3,ob4,ob5,data4=self.delete2dindex(ob1,ob2,ob3,ob4,ob5,data4,index)
            
        
        innercellsY=np.array(np.dstack([ib1,ib2,ib3,ib4,ib5]))
        innercellsY=innercellsY.reshape(ib1.shape[0],5)
        
        outercellsY=np.array(np.dstack([ob1,ob2,ob3,ob4,ob5]))
        outercellsY=outercellsY.reshape(ob1.shape[0],5)
        
        #X-Y Plane boundaries
        count = self.z.shape[0]
        k = np.arange(count)
        
        
        ib2=np.asarray([i for i in k if 0 in z[i] and boundaryx not in x[i] and boundaryy not in y[i]])
        ib1=np.ones(np.size(ib2))*4
        ib3=np.asarray([i+b for i in ib2])
        ib4=np.asarray([i+b+a for i in ib2])
        ib5=np.asarray([i+a for i in ib2])
        data5 = self.image[:,:,0].copy()
        data5[:,:] = 11        
        data5[self.image[:,:,0]==0] = 5
        data5 = data5.flatten()
        index=np.where(data5==11)
        ib1,ib2,ib3,ib4,ib5,data5=self.delete2dindex(ib1,ib2,ib3,ib4,ib5,data5,index)
        #data5=np.ones(ob1.shape[0])*5

        
        ob2=np.asarray([i for i in k if boundaryz in z[i] and boundaryx not in x[i] and boundaryy not in y[i]])
        ob1=np.ones(np.size(ob2))*4
        ob3=np.asarray([i+b for i in ob2])
        ob4=np.asarray([i+b+a for i in ob2])
        ob5=np.asarray([i+a for i in ob2])
        data6 = self.image[:,:,-1].copy()
        data6[:,:] = 12        
        data6[self.image[:,:,-1]==0] = 6
        data6 = data6.flatten()
        #data6=np.ones(ob1.shape[0])*6     
        index=np.where(data6==12)
        ob1,ob2,ob3,ob4,ob5,data6=self.delete2dindex(ob1,ob2,ob3,ob4,ob5,data6,index)
        innercellsZ=np.array(np.dstack([ib1,ib2,ib3,ib4,ib5]))
        innercellsZ=innercellsZ.reshape(ib1.shape[0],5)
        
        outercellsZ=np.array(np.dstack([ob1,ob2,ob3,ob4,ob5]))
        outercellsZ=outercellsZ.reshape(ob1.shape[0],5)
        
        cells=np.vstack([innercellsX,outercellsX,innercellsY,outercellsY,innercellsZ,outercellsZ])
        cellcount=cells.shape[0]
        cells=cells.flatten()
        
                
        quad_type=tvtk.Quad().cell_type
        offset=np.arange(0,cellcount,5)
        
        cell_types=np.ones(cellcount)*quad_type
        
        data1=data1.reshape([np.shape(data1)[0],1])
        data2=data2.reshape([np.shape(data2)[0],1])
        data3=data3.reshape([np.shape(data3)[0],1])
        data4=data4.reshape([np.shape(data4)[0],1])
        data5=data5.reshape([np.shape(data5)[0],1])
        data6=data6.reshape([np.shape(data6)[0],1])
        
        cell_data=np.vstack([data1,data2,data3,data4,data5,data6])
        
        cell_data=cell_data.flatten()
        
        return [cell_types, offset, cellcount, cells, cell_data]
    


    def writeUnstructuredGrid(self,points, cell_types, offset, cell_array, cellData):
        
        ug = tvtk.UnstructuredGrid(points=points)
        ug.set_cells(cell_types, offset, cell_array)
        ug.cell_data.scalars=cellData
        ug.cell_data.scalars.name='MaterialID'
        
        import PythonFCST as fcst
        ob=fcst.mesh.distanceSD(self.image,voxelsize=self.scale,material=0)
        ob.calcPSD(20)
        field_data = ob.imPSD.flatten()
        field_data = np.delete(field_data,np.where(field_data==0))
        field_data = np.append(field_data,self.cellData2)
        ug.cell_data.add_array(field_data)
        ug.cell_data.get_array(1).name='KnudsenRadius'
        ug.cell_data.update()        
        
        writer=tvtk.UnstructuredGridWriter(file_name=self.filename)
        try:
            writer.set_input_data(ug)
        except:
            writer.set_input(ug)
        writer.update()
        writer.write()
    
    def write(self):
        print self._indent, "="*50
        
        points=self.generatePoints()
        [cell_types1, offset1, cellcount1, cells1, cellData1] = self.generate3dCells(points)
        [cell_types2, offset2, cellcount2, cells2, cellData2] = self.generate2dCells(points)
        print self._indent, "="*50
        print self._indent, "=  Points and cells Created"
        cell_types=np.append(cell_types1, cell_types2)
        
        offset=np.append(offset1, offset2)
        
        cellcount= cellcount1 + cellcount2

        cells = np.append(cells1,cells2)
        
        cells = cells.flatten()
        
        cellData=np.append(cellData1, cellData2)
        self.cellData2= cellData2
        cellData=cellData.flatten()
        
        cell_array = tvtk.CellArray()
        cell_array.set_cells(cellcount, cells)
        
        points=np.array(np.dstack([self.x*self.scale[0],self.y*self.scale[1],self.z*self.scale[2]]),'f')
        points=points.reshape(self.x.shape[0],3)
        
        
        print self._indent, "=  Writing a VTK file"
        ug=self.writeUnstructuredGrid(points, cell_types, offset, cell_array, cellData)
        print self._indent, "= The file has been successfully written!"
        print self._indent, "-"*50
        return ug

if __name__ == "__main__":
    image = np.zeros([5,5,2])
    image[0,:,:]=100
    print image.dtype
    o=GridGeneratorLarge(image,'test.vtk')
    o.write()    
