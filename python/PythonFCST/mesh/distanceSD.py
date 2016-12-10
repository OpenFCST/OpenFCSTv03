#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mayank Sabharwal
# License: TBD
# Copyright (c) 2012, TBD


import scipy.ndimage as ndimage
import scipy
import numpy as np
import pylab as pl
import scipy.io
from itertools import cycle


class distanceSD:
    r"""
    Post processor for calculation of pore/particle size distribution based on distance transform method

    image, label='Pore', voxelsize=[1, 1, 1], voxelunit='nm',material=0

    Parameters
    ----------
    image : ndarray
        Input array for the 3D catalyst structure
    label : string
        Specify the material or type of distribution
        Eg: For Pore Size Distribution label = 'Pore'
    voxelsize : tuple
        Pixel size in X,Y,Z direction
        To be specified as [2,2,2] for a pixel size of 2 in each direction
    voxelunit : string
        Unit of length for the voxelsize
    material : int
        Material ID (in the input file) of the component for which the distribution is required


    .. warning:: This documentation is not yet finished.

    .. todo:: Finish documentation

    """
    def __init__(self, image, label='Pore', voxelsize=[1, 1, 1], voxelunit='nm',material=0,**kwargs):

        
        indent='\t'
        print indent, "="*50
        print indent, "Calculatiing Fitted Sphere Size Distribution "
        print indent, "*"*50
        print indent, "= Processing dataset "
        print indent, "-"*50
        self._label=label

        self._spacing=voxelsize

        self._imSolidPore = 0
        self.imDistanceTransform = 0
        self.imPSD = 0
        self._indent=indent
        im = np.asarray(image)
        self._format_color = cycle(['b','g','m','r','y','o'])
        self.scale=np.shape(im)[0]*np.shape(im)[1]
        self._rangeX=scipy.linspace(0,im.shape[0]*voxelsize[0],im.shape[0],endpoint=False)
        self._rangeY=scipy.linspace(0,im.shape[1]*voxelsize[1],im.shape[1],endpoint=False)
        self._rangeZ=scipy.linspace(0,im.shape[2]*voxelsize[2],im.shape[2],endpoint=False)
        self._voxelunit=voxelunit
        self.image=im
        self._imSolidPore=(self.image==material).copy()


    def _getImage(self):
        """
        Private get method: getImage
        """
        return self._image #: Normalized image data

    def _setImage(self, myimage):
        """
        Private set method

        @param myimage: New image
        @type myimage: ndarray
        @return: returns a ndarray
        """
        self._image = myimage.copy() #: Data class of the performance model data
        image = property(fget=_getImage, fset=_setImage) #: Property access  to the data

    def _getImSolidPore(self):
        """
        Private get method: getImage
        """
        return self._imSolidPore #: Normalized image data

    def _setImSolidPore(self, myimage):
        """
        Private set method

        @param myimage: New image
        @type myimage: ndarray
        @return: returns a ndarray
        """
        self._imSolidPore = myimage.copy() #: Data class of the performance model data
        imSolidPore = property(fget=_getImSolidPore, fset=_setImSolidPore) #: Property access  to the data


    def makeDistanceTransform(self,method='edt'):
        """Create Distance transformed image
        Methods for the distance transform:
            - bf    Brute force
            - cdt   Checkerboard
            - edt   Euclidean (default)
        """
        if (method=='edt'):
            self.imDistanceTransform = ndimage.distance_transform_edt(self._imSolidPore,sampling=self._spacing)
        elif (method== 'bf'):
            self.imDistanceTransform = ndimage.distance_transform_bf(self._imSolidPore,sampling=self._spacing)
        elif (method== 'cdt'):
            self.imDistanceTransform = ndimage.distance_transform_cdt(self._imSolidPore,sampling=self._spacing)

    def makePoreSpaceAlg1(self,radius=10.0):
        """
        Finds all the pores where a sphere of radius R can fit into.
        """
        tmp=scipy.ndimage.distance_transform_edt(1-(self.imDistanceTransform>=radius),sampling=self._spacing)
        return ((tmp<=radius)*self._imSolidPore).copy()

    def makePSDAlg1(self,nsteps=10):
        """
        Calculate the comulative pore size distribution. This follows Jeff Gosticks approach described in
        http://10.2.4.201/trac/AFCCNew/wiki/DPM/DPM_IMP/PSD
        """

        print self._indent, "="*50
        print self._indent, "= Calculate PSD"
        span=scipy.linspace(np.min(self._spacing),self.imDistanceTransform.max(),nsteps)
        
        self.imPSD=self.imDistanceTransform.copy()


        for ii in span:
            print self._indent, "=  - Pore radius", ii
            tmp=self.makePoreSpaceAlg1(radius=ii)
            self.imPSD[tmp>0]=ii
      
        return span

    def calcPSD(self, nsteps=10):

        self.makeDistanceTransform()
        self.span=self.makePSDAlg1(nsteps=nsteps)
        psdtmp=[]
        for ii in self.span:
            psdtmp.append(np.average(self.imPSD==ii))
        self.psd=scipy.array(psdtmp)
        self.cpsd = np.cumsum(self.psd)



    def plot_distn(self,nsteps=10, label='Sample1',ax=None,datatype='cpsd', scale='default'):
        self._basename=label
        colors = self._format_color
        if ax==None:
            figure = pl.figure()
            ax = figure.add_subplot(111)
        color=next(colors)
        if datatype == 'cpsd' :
            ax.plot(self.span, self.cpsd,
                        color+'o--',
                        linestyle='--',
                        label=label + ' CPSD'
                        )

        elif datatype =='psd':

            ax.plot(self.span, self.psd,
                        color+'o--',
                        linestyle='--',
                        label=label + ' PSD'
                        )
        pl.xlabel(self._label +' Radius [nm]')
        pl.ylabel(self._label + ' Volume Fraction')
        pl.title('Fiited Sphere Size Distribution')

        pl.legend(loc=1)

        return ax

    def writeCSV(self, filename, datatype='psd'):
        if datatype == 'psd':
            scipy.savetxt(filename,(self.span, self.psd),delimiter=',')
        elif datatype == 'cpsd':
            scipy.savetxt(filename,(self.span, self.cpsd),delimiter=',')
        else:
            raise ValueError