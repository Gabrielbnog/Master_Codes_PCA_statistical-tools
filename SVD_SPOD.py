import scipy
import scipy.linalg
import scipy.signal
from matplotlib import pyplot as plt
from pylab import genfromtxt
import os
from termcolor import colored

#~~~
nsnap = 1537
dt = 0.0005
skip = 2

#~~~
Nf_list = [0,8,16,32]

#~~~
norm_list = ['ru','rv']

#~~~
#mode = 'svd'
#mode = 'full'
mode = 'full'

#~~~
#filter_list = ['gaussian']
filter_list = ['boxcar','gaussian']

#~~~
matrix_path = './SPOD_matrix/'

#~~~
path = './'

#################################################################################################

plt.close('all')

os.system('date')

print ' '

for norm in norm_list:

  print colored ('===================================================', 'cyan' )

  for filter_name in filter_list:
    for Nf in Nf_list:

      print ' '

      char_nf = '%03d' % Nf
      
      output_eigen_path = path + '/output_eigen/' + norm + '_norm/'
      os.system('mkdir -p ' + output_eigen_path)

      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #
      #  
 
      print '    Correlation norm for POD : ', colored ( norm, 'blue' )
      print '  Spectral POD filter family : ', colored ( filter_name, 'red' )
      print '  Spectral POD filter length : ', colored ( char_nf, 'green' )

      full_path = path + matrix_path + '/' + norm + '_norm/'
      filename = 'matrix_' + norm + '_Nf_' + str(char_nf) + '_filter_' + filter_name + '.dat'
      matrix = scipy.zeros((nsnap,nsnap),dtype=float)
      matrix[:,:] = genfromtxt(full_path + filename)

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #
      #  SVD

      print '  SVD'

      U,S,V = scipy.linalg.svd(matrix, full_matrices=1, compute_uv=1)
      eigen_svd = S
      S = scipy.diag(S)

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #
      #  Export temporal modes

      if mode == 'full':

        print '  Export figures'
        
        fig_path = path + 'output_figure/SVD/' + norm + '_Nf_' + char_nf + '_filter_' + filter_name + '/'        

        os.system('mkdir -p ' + fig_path)
        os.system('mkdir -p ' + fig_path + '/phase_diagram')
        os.system('mkdir -p ' + fig_path + '/chronos')
      
        plt.plot(eigen_svd/sum(eigen_svd[:]))
        plt.yscale('log')
        plt.xlim([0,50])
        plt.ylim([4e-3,0.05])
        label = norm + '_eigen_050_svd_Nf_' + str(char_nf) + '_log'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        plt.plot(eigen_svd/sum(eigen_svd[:]))
        plt.yscale('log')
        plt.xlim([0,250])
        plt.ylim([4e-3,0.05])
        label = norm + '_eigen_250_svd_Nf_' + str(char_nf) + '_log'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        plt.plot(eigen_svd/sum(eigen_svd[:]))
        plt.xlim([0,50])
        plt.ylim([0,1.1*max(eigen_svd)/sum(eigen_svd[:])])
        label = norm + '_eigen_050_svd_Nf_' + str(char_nf) + '_lin'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        plt.plot(eigen_svd/sum(eigen_svd[:]))
        plt.xlim([0,250])
        plt.ylim([0,1.1*max(eigen_svd)/sum(eigen_svd[:])])
        label = norm + '_eigen_250_svd_Nf_' + str(char_nf) + '_lin'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        sum_modes = scipy.zeros(nsnap)
        for i in range(0,nsnap):
          sum_modes[i] = (eigen_svd[i]+sum(eigen_svd[:i]))/sum(eigen_svd[:])
          
        plt.plot(1.0-sum_modes,'+',color='blue')
        plt.plot(1.0-sum_modes,color='blue')
        plt.plot(scipy.ones(nsnap)*0.01,'--',color='red',linewidth=0.75)  
        plt.plot(scipy.ones(nsnap)*0.05,'--',color='orange',linewidth=0.75)
        plt.yscale('log')
        plt.xlim([0,50])
        plt.ylim([9e-3,1.0])
        label = norm + '_sum_eigen_050_svd_Nf_' + str(char_nf) + '_log'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()
         
        plt.plot(1.0-sum_modes,'+',color='blue')
        plt.plot(1.0-sum_modes,color='blue')
        plt.plot(scipy.ones(nsnap)*0.01,'--',color='red',linewidth=0.75)  
        plt.plot(scipy.ones(nsnap)*0.05,'--',color='orange',linewidth=0.75)
        plt.yscale('log')
        plt.xlim([0,250])
        plt.ylim([9e-3,1.0])
        label = norm + '_sum_eigen_250_svd_Nf_' + str(char_nf) + '_log'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        plt.plot(sum_modes,'+',color='blue')
        plt.plot(sum_modes,color='blue')
        plt.plot(scipy.ones(nsnap)*0.99,'--',color='red',linewidth=0.75)  
        plt.plot(scipy.ones(nsnap)*0.95,'--',color='orange',linewidth=0.75)  
        plt.xlim([0,50])
        plt.ylim([0,1.0])
        label = norm + '_sum_eigen_050_svd_Nf_' + str(char_nf) + '_lin'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        plt.plot(sum_modes,'+',color='blue')
        plt.plot(sum_modes,color='blue')
        plt.plot(scipy.ones(nsnap)*0.99,'--',color='red',linewidth=0.75)  
        plt.plot(scipy.ones(nsnap)*0.95,'--',color='orange',linewidth=0.75)  
        plt.xlim([0,250])
        plt.ylim([0,1.0])
        label = norm + '_sum_eigen_250_svd_Nf_' + str(char_nf) + '_lin'
        plt.savefig(fig_path + label + '.pdf')
        plt.close()

        time = scipy.zeros(nsnap)
        for k in range(0,nsnap):
          time[k] = float(k)*float(skip)*dt

        for k in range(0,60,2):
          plt.plot(time[:],U[:,k])
          plt.plot(time[:],U[:,k+1])
          #plt.plot(U[:,k]*scipy.sqrt(eigen_svd[k]))
#          plt.ylim([-1,1])
          k1 = k+1
          char_k = '%02d' % k1
          label = norm + '_chronos_' + char_k
          #plt.savefig(fig_path + '/chronos/' + label + '.pdf')
          plt.savefig(fig_path + '/chronos/' + label + '.png',dpi=600)
          plt.close()

        for k in range(0,59,2):
          k1 = k+1
          k2 = k+2
          char_1 = '%02d' % k1
          char_2 = '%02d' % k2
          plt.plot(U[:,k],U[:,k+1])
          plt.axis('equal')
          label = norm + '_phase_diagram_' + char_1 + '_' + char_2
          #plt.savefig(fig_path + '/phase_diagram/' + label + '.pdf')
          plt.savefig(fig_path + '/phase_diagram/' + label + '.png',dpi=600)
          plt.close()

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #
      #  Export matrix

      scipy.savetxt( path + '/output_eigen/' + norm + '_norm/eigenvalues_' + norm + '_Nf_' + char_nf + '_filter_' + filter_name + '.dat', eigen_svd )
      scipy.savetxt( path + '/output_eigen/' + norm + '_norm/eigenvector_' + norm + '_Nf_' + char_nf + '_filter_' + filter_name + '.dat', U, fmt='%25.18e' )

print ' '
print 'End'
