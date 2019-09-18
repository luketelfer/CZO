# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:37:03 2016

@author: miguel

modified pfb writer, for soil,vege, and parflow input files that are in pfb format
by Chao Chen, Oct 29, 2018
"""
import numpy as np
import struct
from numpy import genfromtxt
import multiprocessing as mp
import time





###################################
###################################
# function, write the data into pfb format
###################################


###################################
def pfb_writer_multilayer(x1, y1, z1, nx, ny, nz, dx, dy, dz, V1, ix, iy, iz, nnx, nny, nnz, rx, ry, rz, ns, fn):
    with open(fn, 'wb') as fid:
        X = struct.pack('>d', x1)
        fid.write(X)
        Y = struct.pack('>d', y1)
        fid.write(Y)
        Z = struct.pack('>d', z1)
        fid.write(Z)

        NX = struct.pack('>i', nx)
        fid.write(NX)
        NY = struct.pack('>i', ny)
        fid.write(NY)
        NZ = struct.pack('>i', nz)
        fid.write(NZ)

        DX = struct.pack('>d', dx)
        fid.write(DX)
        DY = struct.pack('>d', dy)
        fid.write(DY)
        DZ = struct.pack('>d', dz)
        fid.write(DZ)

        NS = struct.pack('>i', ns)
        fid.write(NS)

        for i0 in range(ns):
            iX = struct.pack('>i', ix)
            fid.write(iX)
            iY = struct.pack('>i', iy)
            fid.write(iY)
            iZ = struct.pack('>i', iz)
            fid.write(iZ)
            nnX = struct.pack('>i', nnx)
            fid.write(nnX)
            nnY = struct.pack('>i', nny)
            fid.write(nnY)
            nnZ = struct.pack('>i', nnz)
            fid.write(nnZ)
            rX = struct.pack('>i', rx)
            fid.write(rX)
            rY = struct.pack('>i', ry)
            fid.write(rY)
            rZ = struct.pack('>i', rz)
            fid.write(rZ)

            for k in range(iz, iz + nnz):
                # print 'k = %i' %k
                for i in range(iy, iy + nny):
                    # print 'i = %i' %i
                    for j in range(ix, ix + nnx):
                        # print 'j = %i' %j
                        VAR = struct.pack('>d',V1[i,j,k])
                        fid.write(VAR)


def pfb_writer_onelayer(x1,y1,z1,nx,ny,nz,dx,dy,dz,V1,ix,iy,iz,nnx,nny,nnz,rx,ry,rz,ns,fn):

   with open(fn,'wb') as fid:
        X = struct.pack('>d',x1)
        fid.write(X)
        Y = struct.pack('>d',y1)
        fid.write(Y)
        Z = struct.pack('>d',z1)
        fid.write(Z)
#        VAR = np.empty(shape=(np.int(ny),np.int(nx),np.int(nz)))
        
        NX = struct.pack('>i',nx)
        fid.write(NX)
        NY = struct.pack('>i',ny)
        fid.write(NY)
        NZ = struct.pack('>i',nz)
        fid.write(NZ)

        DX = struct.pack('>d',dx)
        fid.write(DX)
        DY = struct.pack('>d',dy)
        fid.write(DY)
        DZ = struct.pack('>d',dz)
        fid.write(DZ)

        NS = struct.pack('>i',ns)
        fid.write(NS)
#        print V1
        for i0 in range(ns):
            iX = struct.pack('>i',ix)
            fid.write(iX)
            iY = struct.pack('>i',iy)
            fid.write(iY)
            iZ = struct.pack('>i',iz)
            fid.write(iZ)
            nnX = struct.pack('>i',nnx)
            fid.write(nnX)
            nnY = struct.pack('>i',nny)
            fid.write(nnY)
            nnZ = struct.pack('>i',nnz)
            fid.write(nnZ)
            rX = struct.pack('>i',rx)
            fid.write(rX)
            rY = struct.pack('>i',ry)
            fid.write(rY)
            rZ = struct.pack('>i',rz)
            fid.write(rZ)
                     
            for k in range(iz,iz+nnz):
                for i in range(iy,iy+nny):
                    for j in range(ix,ix+nnx):
#                        print V1[i,j]
                        VAR = struct.pack('>d',V1[i,j])
#                        print VAR
                        #fid.write(struct.pack('>d',VAR[i,j,k]))
                        fid.write(VAR)
       
