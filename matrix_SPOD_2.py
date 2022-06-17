import scipy
import scipy.linalg
import scipy.signal
from matplotlib import pyplot as plt
from pylab import genfromtxt
import os
from termcolor import colored
import glob
import time

#~~~ Number of snapshots (.cfl file)
nsnap = 1537

#~~~ Half size of the filter
half_filter_list = [0,8,16,32]

#~~~ Periodicity of the correlation matrix
periodicity = 1

#~~~ path to the matrix of each zone
path = './corr_matrix/'

#~~~ Name of the files containing the correlation matrix
#pattern = '/corr_matrix_????_'

#~~~ Name of the SPOD filter
filter_list = ['boxcar','gaussian']
#filter_list = ['boxcar']
#filter_list = ['gaussian']

#~~~ POD norm
#norm_list = ['tke','ru','rv','rw']
norm_list = ['ru','rv']

#################################################################################################

os.system('date')

for norm in norm_list:

  print ' '
  print ' ========================================================================== '  
  print ' '  
  
  if norm == 'prs':
    print colored( '  + Norma de Pressao', 'blue' )
    variable_list = 'p'
  if norm == 'tke':
    print colored( '  + Norma de Energia Cinetica', 'yellow' )
    variable_list = 'k'
  if norm == 'ru':
    print colored( '  + Norma de Momento X', 'red' )
    variable_list = 'ru'
  if norm == 'rv':
    print colored( '  + Norma de Momento Y', 'green' )
    variable_list = 'rv'      
  if norm == 'rw':
    print colored( '  + Norma de Momento Z', 'cyan' )
    variable_list = 'rw'

  output_figure_path = './output_figure/corr_matrix/' + norm + '_norm/'
  output_matrix_path = './SPOD_matrix/' + norm + '_norm/'
  os.system('mkdir -p ' + output_figure_path)
  os.system('mkdir -p ' + output_matrix_path)

  for filter_name in filter_list:

    for Nf in half_filter_list:

      #    
      start = time.time()

      print ' '
      
      char_nf = '%03d' % Nf

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    #  Matrix input

      matrix_file_unfiltered = output_matrix_path + 'matrix_' + norm + '_Nf_000.dat'
      if os.path.isfile(matrix_file_unfiltered):

        matrix = genfromtxt(matrix_file_unfiltered)

      else:

        matrix = genfromtxt(path + 'corr_matrix_' + variable_list + '.out')

#        matrix = scipy.zeros((nsnap,nsnap))
#        matrix_input = scipy.zeros((nsnap,nsnap))
#        
#        #~~~
#        i = 0
#        for filename in list_of_files:
#          
#          i = i + 1
#          print i, filename
#          matrix_input[:,:] = genfromtxt(filename)
#          
#          #checking for NaN
#          if (matrix[0,0] != matrix[0,0]):
#            print 'skipping' + filename
#          else:
#            matrix = matrix + matrix_input[:,:]        

#        print ' '

      #~~~
      plt.figure()
      plt.contourf(matrix)
      plt.xlim([0,nsnap])
      plt.ylim([0,nsnap])
      plt.colorbar()
      plt.savefig(output_figure_path + norm + '_full_matrix_spod_Nf_000' + '_filter_' + filter_name + '.pdf')
      plt.savefig(output_figure_path + norm + '_full_matrix_spod_Nf_000' + '_filter_' + filter_name + '.png',dpi=600)   
      plt.close()

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    #  Spectral POD

      #~~~ Snapshot POD
      if Nf == 0:
        print '  Spectral POD filter :', colored( Nf, 'red')
        print '    -> No filtering employed ...'

      #~~~ Evitar bizonhice
      if 2*Nf+1 > nsnap:
        print 'O filtro eh maior que a matriz!'
        break
     
      #~~~ Filter set-up
      if Nf > 0:
        print '  Spectral POD filter :', colored( Nf, 'red')
        print '  Total filter width = (2N+1) =', 2*Nf+1
        if filter_name == 'gaussian':    
          print '    -> Filter : ', colored( 'Gaussian', 'cyan' )
          sigma = Nf/4.0
          g = scipy.zeros(2*Nf+1)
          k = scipy.zeros(2*Nf+1)
          for i in range (0,2*Nf+1):
            k[i] = i-Nf
            g[i] = scipy.exp(-0.5*(k[i]/sigma)**2)
        if filter_name == 'boxcar':
          print '    -> Filter : ', colored( 'boxcar', 'cyan' )
          g = scipy.signal.boxcar(2*Nf+1,sym=True)
        #if filter_name == 'kaiser14':
        #  print '    -> Filter : Kaiser_14'
        #  g = scipy.signal.kaiser(2*Nf+1,14.0,sym=True)
     
      #~~~ Filtering
      if Nf > 0:
      
        if periodicity == 0:
        
          print '  No periodicity.. '
          
          matrix_bkp = matrix

          for i in range(0,nsnap):
            for j in range(i,nsnap):
              gaux = scipy.zeros(2*Nf+1)
              dummy = 0.0
              cont = 0
              for l in range(0,2*Nf+1):
                k = l-Nf
                if (i+k >= 0 and i+k < nsnap):
                  if (j+k >= 0 and j+k < nsnap):
                    cont = cont + 1
                    dummy = dummy + matrix_bkp[i+k,j+k]*g[l]
                    gaux[l] = g[l]
              if cont == 0:
                print 'FUDEU!!'
              if sum(scipy.absolute(gaux[:])) < 0.000001:
                print 'FUDEU!!'
              matrix[i,j] = dummy/sum(gaux[:])
            for j in range(i,nsnap):
              matrix[j,i] = matrix[i,j]

        if periodicity == 1:
        
          print '  Using periodicity.. Stacking the matrix itself!'
          dummy1 = scipy.hstack((matrix,matrix,matrix))
          stack_matrix = scipy.vstack((dummy1,dummy1,dummy1))
        
          for i in range(0,nsnap):
            for j in range(i,nsnap):
              dummy = 0.0
              for l in range(0,2*Nf+1):
                k = l-Nf
                dummy = dummy + stack_matrix[i+k+nsnap,j+k+nsnap]*g[l]
              matrix[i,j] = dummy/sum(g[:])
            for j in range(i,nsnap):
              matrix[j,i] = matrix[i,j]

        plt.figure()
        plt.contourf(matrix)
        plt.xlim([0,nsnap])
        plt.ylim([0,nsnap])
        plt.colorbar()
        label = norm + '_full_matrix_spod_Nf_' + str(char_nf) + '_filter_' + filter_name
        plt.savefig(output_figure_path + label + '.pdf')
        plt.savefig(output_figure_path + label + '.png',dpi=600)
        plt.close()
                  
      plt.close('all')

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #
      #  Export matrix
      
      if Nf == 0:
        scipy.savetxt(output_matrix_path + '/matrix_' + norm + '_Nf_000.dat', matrix, fmt='%25.18e')
      scipy.savetxt(output_matrix_path + '/matrix_' + norm + '_Nf_' + str(char_nf) + '_filter_' + filter_name + '.dat', matrix, fmt='%25.18e')
      
      #
      end = time.time()
      print '  Elapsed time : ', (end - start), 'seconds'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print ' '
print 'End'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
