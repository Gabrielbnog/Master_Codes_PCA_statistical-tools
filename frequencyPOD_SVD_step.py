
#precisa checar nesse codigo como eh feita a leitura dos arquivos intermediarios que saem do xPOD
#precisa ser feitas algumas manipulacoes, eu acho. Mas se precisar, eh soh usando funcao do Linux
#precisa ter cuidado com o overlap dos pontos de malha para nao cagar o SVD

import numpy as np
from pylab import genfromtxt
import matplotlib.pyplot as plt
import scipy
from scipy import linalg
from scipy.io import FortranFile

import sys

np.set_printoptions(edgeitems=10,linewidth=132,suppress=True)

#~~~~~~~~~~~~~~~~~~~~~~~

join_zones = True
single_zone = False

if join_zones == True:

  nzones = 5
  header = np.zeros((nzones+1,5),dtype=np.int32)
  header[1:nzones+1,:] = genfromtxt('header.dat')

  header[0,1] = np.sum(header[1:,1])
  header[0,2] = header[1,2]
  header[0,3] = header[1,3]
  header[0,4] = header[1,4]

  imin = np.zeros((nzones+1),dtype=np.int32)
  imax = np.zeros((nzones+1),dtype=np.int32)

  imax[0] = 0
  for m in range(1,nzones+1):

    imin[m] = imax[m-1] + 1
    imax[m] = imin[m] + header[m,1] - 1
    
    print imin[m], imax[m]
    
  imin[0] = 1
  imax[0] = imax[-1]

  print imin[0], imax[0]
  
  nx1 = header[0,1]
  nx2 = header[0,2]
  nx3 = header[0,3]
  ngrid = nx1*nx2*nx3
  
  nbins = header[0,4]
  
  print ngrid, nbins, '\n'
  
  n = 0
  grid_file = open('dummy_grid.dat','w')
  grid_file.write('FILETYPE=GRID' + '\n')
  grid_file.write('VARIABLES="X","Y","Z"' + '\n')
  for m in range(1,nzones+1):

    char_zone = '%02d' % m
    
    print ' Working in zone ' + char_zone

    nx1z = header[m,1]
    nx2z = header[m,2]
    nx3z = header[m,3]
    ngridz = nx1z*nx2z*nx3z

    print nx1z, nx2z, nx3z

    #~~~~~~~~~~~~~~~~~~~~~~~

    print ' Leitura da malha'
    
    xcoord = np.zeros(ngrid)
    ycoord = np.zeros(ngrid)
    zcoord = np.zeros(ngrid)        
    
    # PRECISA ter esses arquivos abaixo!
    xcoord = genfromtxt('./grid/grid_xcoord_z' + char_zone + '.dat')
    ycoord = genfromtxt('./grid/grid_ycoord_z' + char_zone + '.dat')
    zcoord = genfromtxt('./grid/grid_zcoord_z' + char_zone + '.dat')

    grid_file.write('ZONE T="z' + char_zone + '",I=' + str(nx1z) + ',J=' + str(nx2z) + ',K=' + str(nx3z) + ',F=BLOCK' + '\n')

    print 'x'
    idx = 0
    for k in range(0,nx3z):
      for j in range(0,nx2z):
        for i in range(0,nx1z):
          grid_file.write(str(xcoord[idx]) + '\n')
          idx = idx + 1

    print 'y'
    idx = 0
    for k in range(0,nx3z):
      for j in range(0,nx2z):
        for i in range(0,nx1z):
          grid_file.write(str(ycoord[idx]) + '\n')
          idx = idx + 1

    print 'z'
    idx = 0
    for k in range(0,nx3z):
      for j in range(0,nx2z):
        for i in range(0,nx1z):
          grid_file.write(str(zcoord[idx]) + '\n')
          idx = idx + 1

    print idx

#    #~~~~~~~~~~~~~~~~~~~~~~~

#    print ' Leitura dos dados - zonal'
#    
#    aux = np.zeros((ngridz,2*nbins))
#    aux = genfromtxt('./fourier_matrix_' + char_zone + '.dat')

#    #~~~~~~~~~~~~~~~~~~~~~~~
#    print ' Convertendo para complexo ... \n'
#    print n
#    for i in range(0,ngridz):
#      if (i != 0):
#        if (i%500000 == 0):
#          print i
#      for k in range(0,nbins):
#        data[n,k] = aux[i,2*k] + 1j*aux[i,2*k+1]
#      n = n + 1
#      
#    print n

  grid_file.close()
  
  #~~~~~~~~~~~~~~~~~~~~~~~

  print ' Alocando memoria ...'
  
  aux = np.zeros((ngrid,2*nbins))
  
  data = np.zeros((ngrid,nbins),dtype=np.complex128)
  
  #~~~~~~~~~~~~~~~~~~~~~~~

  print ' Leitura dos dados - global'

  # PRECISA ter esse arquivo!!
  aux = genfromtxt('./fourier_matrix.dat')

  #~~~~~~~~~~~~~~~~~~~~~~~
  print ' Convertendo para complexo ... \n'
  
  for i in range(0,ngrid):
    if (i != 0):
      if (i%512000 == 0):
        print '   ', i
    for k in range(0,nbins):
      data[i,k] = aux[i,2*k] + 1j*aux[i,2*k+1]

  #~~~~~~~~~~~~~~~~~~~~~~~
  print ' SVD dos dados ...'
  Ux,Sx,Vxt = scipy.linalg.svd(data, full_matrices=0, compute_uv=1)
  Vx = np.transpose(np.conjugate(Vxt))

  print Ux.shape
  print Sx.shape
  print Vx.shape

  #~~~~~~~~~~~~~~~~~~~~~~~
  print ' Exporting ...'
  
  np.savetxt('lambda.dat',Sx)

  mode_file = open('spatial_16_modes_SPOD__IMAG__.dat','w')

  mode_file.write('FILETYPE=SOLUTION' + '\n')
  for n in range(0,16):
    char_mode = '%02d' % (n+1)
    mode_file.write('VARIABLES="mode' + char_mode + '"'  + '\n')

  idx_f = 0
  for m in range(1,nzones+1):

    char_zone = '%02d' % m
    
    print ' Exporting zone ' + char_zone
  
    nx1z = header[m,1]
    nx2z = header[m,2]
    nx3z = header[m,3]
  
    mode_file.write('ZONE T="z' + char_zone + '",I=' + str(nx1z) + ',J=' + str(nx2z) + ',K=' + str(nx3z) + ',F=BLOCK' + '\n')
   
    for n in range(0,16):
      char_mode = '%02d' % (n+1)
      idx = idx_f
      print 'mode ' + char_mode
      for k in range(0,nx3z):
        for j in range(0,nx2z):
          for i in range(0,nx1z):
            #mode_file.write(str(np.real(Ux[idx,n])) + '\n')
            mode_file.write(str(np.imag(Ux[idx,n])) + '\n')
            idx = idx + 1

    for k in range(0,nx3z):
      for j in range(0,nx2z):
        for i in range(0,nx1z):
          idx_f = idx_f + 1

  mode_file.close()

print 'Fim! (1)'
sys.exit()


#######################################################################3







#######################################################################3
if single_zone == True:
    
  nzones = 5
  header = np.zeros((nzones+1,5),dtype=np.int32)
  header[1:nzones+1,:] = genfromtxt('header.dat')

  for m in range(1,nzones+1):

    char_zone = '%02d' % m
    
    print ' Working in zone ' + char_zone

    nx1 = header[m,1]
    nx2 = header[m,2]
    nx3 = header[m,3]
    nbins = header[m,4]
    
    ngrid = nx1*nx2*nx3

    #~~~~~~~~~~~~~~~~~~~~~~~

    print ' Leitura da malha'
    
    xcoord = np.zeros(ngrid)
    ycoord = np.zeros(ngrid)
    zcoord = np.zeros(ngrid)        
    
    xcoord = genfromtxt('./grid/grid_xcoord_z' + char_zone + '.dat')
    ycoord = genfromtxt('./grid/grid_ycoord_z' + char_zone + '.dat')
    zcoord = genfromtxt('./grid/grid_zcoord_z' + char_zone + '.dat')
    
    mode_file = open('dummy_z' + char_zone + '.dat','w')

    mode_file.write('VARIABLES="X","Y","Z"' + '\n')
    mode_file.write('ZONE T="mode1",I=' + str(nx1) + ',J=' + str(nx2) + ',K=' + str(nx3) + ',F=BLOCK' + '\n')

    print 'x'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(xcoord[idx]) + '\n')
          idx = idx + 1

    print 'y'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(ycoord[idx]) + '\n')
          idx = idx + 1

    print 'z'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(zcoord[idx]) + '\n')
          idx = idx + 1

    mode_file.close()      

    #~~~~~~~~~~~~~~~~~~~~~~~

    print ' Leitura dos dados'

    aux = np.zeros((ngrid,2*nbins))
    aux = genfromtxt('./fourier_matrix_' + char_zone + '.dat')

    #~~~~~~~~~~~~~~~~~~~~~~~
    print ' Convertendo para complexo ...'
    data = np.zeros((ngrid,nbins),dtype=np.complex128)
    for i in range(0,ngrid):
      if (i%100000 == 0):
        print i
      for k in range(0,nbins):
        data[i,k] = aux[i,2*k] + 1j*aux[i,2*k+1]

    #~~~~~~~~~~~~~~~~~~~~~~~
    print ' SVD dos dados ...'
    Ux,Sx,Vxt = scipy.linalg.svd(data, full_matrices=0, compute_uv=1)
    Vx = np.transpose(np.conjugate(Vxt))

    print Ux.shape
    print Sx.shape
    print Vx.shape

    #~~~~~~~~~~~~~~~~~~~~~~~
    print ' Exporting ...'  
    
    np.savetxt('lambda.dat',Sx)

    mode_file = open('spatial_mode_4_modes_SPOD_z' + char_zone + '.dat','w')

    mode_file.write('VARIABLES="X","Y","Z","mode1","mode2","mode3","mode4"' + '\n')
    mode_file.write('ZONE T="mode1",I=' + str(nx1) + ',J=' + str(nx2) + ',K=' + str(nx3) + ',F=BLOCK' + '\n')

    print 'x'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(xcoord[idx]) + '\n')
          idx = idx + 1

    print 'y'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(ycoord[idx]) + '\n')
          idx = idx + 1

    print 'z'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(zcoord[idx]) + '\n')
          idx = idx + 1

  #  mode_file = open('spatial_mode_n01_SPOD.dat','w')

  #  mode_file.write('FILETYPE=SOLUTION' + '\n')
  #  mode_file.write('VARIABLES="mode1"' + '\n')
  #  mode_file.write('ZONE T="mode1",I=' + str(nx1) + ',J=' + str(nx2) + ',K=' + str(nx3) + ',F=BLOCK' + '\n')

    print 'mode01'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(np.real(Ux[idx,0])) + '\n')
          idx = idx + 1
          
    print 'mode02'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(np.real(Ux[idx,1])) + '\n')
          idx = idx + 1
          
    print 'mode03'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(np.real(Ux[idx,2])) + '\n')
          idx = idx + 1
          
    print 'mode04'
    idx = 0
    for k in range(0,nx3):
      for j in range(0,nx2):
        for i in range(0,nx1):
          mode_file.write(str(np.real(Ux[idx,3])) + '\n')
          idx = idx + 1     

    mode_file.close()    







#if join_zones == True:

#  nzones = 5
#  header = np.zeros((nzones+1,5),dtype=np.int32)
#  header[1:nzones+1,:] = genfromtxt('header.dat')

#  header[0,1] = np.sum(header[1:,1])
#  header[0,2] = header[1,2]
#  header[0,3] = header[1,3]
#  header[0,4] = header[1,4]

#  nx1 = header[0,1]
#  nx2 = header[0,2]
#  nx3 = header[0,3]
#  ngrid = nx1*nx2*nx3
#  
#  nbins = header[0,4]
# 
#  aux = np.zeros((ngrid,2*nbins))
#  data = np.zeros((ngrid,nbins),dtype=np.complex128)
#  
#  print ngrid, nbins, '\n'
#  
#  print nx1, nx2, nx3

#  #~~~~~~~~~~~~~~~~~~~~~~~

#  print ' Leitura da malha'
#  
#  xcoord = np.zeros(ngrid)
#  ycoord = np.zeros(ngrid)
#  zcoord = np.zeros(ngrid)        
#  
#  xcoord = genfromtxt('./grid/grid_xcoord.dat')
#  ycoord = genfromtxt('./grid/grid_ycoord.dat')
#  zcoord = genfromtxt('./grid/grid_zcoord.dat')

#  #~~~~~~~~~~~~~~~~~~~~~~~

#  print ' Leitura dos dados'
#  
#  aux = genfromtxt('./fourier_matrix.dat')

#  #~~~~~~~~~~~~~~~~~~~~~~~
#  print ' Convertendo para complexo ... \n'
#  print n
#  for i in range(0,ngrid):
#    if (i != 0):
#      if (i%500000 == 0):
#        print i
#    for k in range(0,nbins):
#      data[i,k] = aux[i,2*k] + 1j*aux[i,2*k+1]

#  #~~~~~~~~~~~~~~~~~~~~~~~
#  print ' SVD dos dados ...'
#  Ux,Sx,Vxt = scipy.linalg.svd(data, full_matrices=0, compute_uv=1)
#  Vx = np.transpose(np.conjugate(Vxt))

#  print Ux.shape
#  print Sx.shape
#  print Vx.shape

#  #~~~~~~~~~~~~~~~~~~~~~~~
#  print ' Exporting ...'
#  
#  np.savetxt('lambda.dat',Sx)

#  mode_file = open('spatial_mode_n01_SPOD.dat','w')

##  mode_file.write('FILETYPE=SOLUTION' + '\n')
##  mode_file.write('VARIABLES="mode1"' + '\n')
##  mode_file.write('ZONE T="mode1",I=' + str(nx1) + ',J=' + str(nx2) + ',K=' + str(nx3) + ',F=BLOCK' + '\n')
#  
#  mode_file.write('VARIABLES="X","Y","Z","mode01"' + '\n')
#  mode_file.write('ZONE T="mode1",I=' + str(nx1) + ',J=' + str(nx2) + ',K=' + str(nx3) + ',F=BLOCK' + '\n')

#  print 'x'
#  idx = 0
#  for k in range(0,nx3):
#    for j in range(0,nx2):
#      for i in range(0,nx1):
#        mode_file.write(str(xcoord[idx]) + '\n')
#        idx = idx + 1

#  print 'y'
#  idx = 0
#  for k in range(0,nx3):
#    for j in range(0,nx2):
#      for i in range(0,nx1):
#        mode_file.write(str(ycoord[idx]) + '\n')
#        idx = idx + 1

#  print 'z'
#  idx = 0
#  for k in range(0,nx3):
#    for j in range(0,nx2):
#      for i in range(0,nx1):
#        mode_file.write(str(zcoord[idx]) + '\n')
#        idx = idx + 1

#  print 'mode01'
#  idx = 0
#  for k in range(0,nx3):
#    for j in range(0,nx2):
#      for i in range(0,nx1):
#        mode_file.write(str(np.real(Ux[idx,0])) + '\n')
#        idx = idx + 1

#  mode_file.close()

##print 'Fim! (1)'
##sys.exit()
