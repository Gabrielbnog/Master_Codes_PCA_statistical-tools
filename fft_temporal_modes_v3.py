
import scipy
import scipy.signal
import numpy.fft
import numpy
from pylab import genfromtxt
from matplotlib import pyplot as plt
import os
import matplotlib.ticker as plticker

nsnap = 1537
dt = 0.0005
skip = 100

Nf_list=['000','008','016','032','128','256','512']
#Nf_list=['000','256']

norm_list = ['ru','rv']

nmodes = 40

filter_list = ['gaussian','boxcar']

desired_frequency = [9.0]

FFT = True
Sigma = False

path = './output_eigen/'
path_to_figures = './output_figure/'

#########################################################################

plt.close('all')

os.system('date')

for filter_name in filter_list:

  for norm in norm_list:

    if norm == 'prs': print 'Norma de Pressao'
    if norm == 'tke': print 'Norma de Energia Cinetica'
    if norm == 'ru' : print 'Norma de Momento U'
    if norm == 'rv' : print 'Norma de Momento V'
    if norm == 'rw' : print 'Norma de Momento W'

    if filter_name == 'gaussian' : print ' Gaussian filter'
    if filter_name == 'boxcar' : print ' Boxcar filter'

    path_output_figures = path_to_figures + '/FFT/' + norm + '_fft_modes_filter_' + filter_name + '/'
    os.system('mkdir -p ' + path_output_figures)

    U = scipy.zeros((nsnap,nsnap,len(Nf_list)))
    S = scipy.zeros((nsnap,len(Nf_list)))
    for k in range(0,len(Nf_list)):
      U[:,:,k] = genfromtxt(path + norm + '_norm/eigenvector_' + norm + '_Nf_' + Nf_list[k] + '_filter_' + filter_name + '.dat')
      S[:,k]   = genfromtxt(path + norm + '_norm/eigenvalues_' + norm + '_Nf_' + Nf_list[k] + '_filter_' + filter_name + '.dat')

      for i in range(0,nsnap):
        U[:,i,k] = U[:,i,k]*scipy.sqrt(S[i,k])
      
    #scipy.savetxt('frequency.dat', xf[:], fmt='%25.18f')

    #########################################################################

    if FFT == True:
      ###################
      #
      #  FFT computation

      freq = scipy.ones((2*len(desired_frequency),2))
      for i in range(len(desired_frequency)):
        k1 = 2*i
        k2 = k1 + 1
        freq[k1,0] = desired_frequency[i]
        freq[k2,0] = freq[k1,0]
        freq[k1,1] = 1e-16
        freq[k2,1] = 1e+5

      N = nsnap/2
      print N
      Pxx = scipy.zeros((N/2+1,len(Nf_list)))
      for i in range(0,nmodes,2):

        ii = i + 1
        char = '%03d' % ii

        window_name = 'hamming'

        for k in range(0,len(Nf_list)):
        
          x = U[:,i,k]
          [Pxx[:,k], freqs] = plt.psd(x, NFFT=N, Fs=20.0*2.0*numpy.pi, Fc=0, detrend='mean',
             window=numpy.hamming(N), noverlap=2*N/3, pad_to=None,
             sides='onesided', scale_by_freq=True, return_line=None)

        plt.close()
        
        fig = plt.figure(1)
        #fig.suptitle('Flyover', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        #fig.subplots_adjust(top=0.85)
        ax.set_title(filter_name + ' - ' + window_name + ' - mode : ' + char, fontsize=13, fontweight='bold',position=(0.5,1.02))  
        
        for k in range(0,len(Nf_list)):    
          name = Nf_list[k]
          plt.plot(freqs[:],Pxx[:,k],marker='.',markersize=3.0,label=str(name))
          for l in range(0,len(desired_frequency)):
            plt.plot(freq[2*l:2*l+2,0],freq[2*l:2*l+2,1],'--',color='black',linewidth=0.75)

  #      ax.set_position([0.14,0.14,0.72,0.72])
  #  #    plt.xlim([0000,1500])
  #  #    plt.ylim([1e-4,5e+4])

  #      plt.xlim([0,60])

        plt.ylim([1e-12,1e-5])    
        if norm == 'prs':
          plt.ylim([1e-15,1e-3])

        plt.xlim([0,40])
  #      plt.xscale('log')
        plt.yscale('log')
  #      plt.grid('on')
        plt.minorticks_on()
        ax.xaxis.tick_top()
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.tick_params(labelbottom='on',labeltop='off')
  #      loc = plticker.MultipleLocator(base=300) # this locator puts ticks at regular intervals
  #      ax.xaxis.set_major_locator(loc)  
  #      ax.tick_params(axis='x', labelsize=12)
  #      ax.tick_params(axis='y', labelsize=12)
  #      plt.ylabel('', fontsize=13, fontweight='bold')
  #      plt.xlabel('Frequency (Hz)', fontsize=13, fontweight='bold')
        plt.legend(loc='upper right',fontsize=10)
  #      #ax.get_xaxis().set_major_formatter(plticker.ScalarFormatter())
  #      #ax.set_xticks([100, 500, 1000, 5000, 10000]) 
  #      #ax.ticklabel_format(axis='x', style='sci')
  #      ax.xaxis.labelpad = 12  #skip entre o label e o grafico
  #      ax.yaxis.labelpad = 12  #skip entre o label e o g
        ax.text(11.0, 0.000025, 'kc ~ 9', fontsize=10,fontweight='bold',bbox={'facecolor':'white', 'alpha':0.33, 'pad':4})        
  #  #    ax.text(275.0, 6000.0, '240Hz', fontsize=10,fontweight='bold',bbox={'facecolor':'white', 'alpha':0.33, 'pad':4})
  #  #    ax.text(500.0, 1000.0, '470Hz', fontsize=10,fontweight='bold',bbox={'facecolor':'white', 'alpha':0.33, 'pad':4})
  #  #    ax.text(975.0, 0.15, '940Hz', fontsize=10,fontweight='bold',bbox={'facecolor':'white', 'alpha':0.33, 'pad':4})
  #  #    ax.text(1950.0, 1000.0, '1950Hz', fontsize=10,fontweight='bold',bbox={'facecolor':'white', 'alpha':0.33, 'pad':4})
  #  #    plt.savefig(path_output_figures + '/fft' + char + '.pdf')
        plt.savefig( path_output_figures + '/fft' + char + '.png',dpi=750)
        plt.close()
        #plt.show()

    #########################################################################
      
    ###################
    #
    #  Eigenvalues

    if Sigma == True:
      path_output_figures = path + 'output_figure/SVD/'
      #os.system('mkdir -p ' + path_output_figures)

      for k in range(0,len(Nf_list)):
        print 'The sum of the eigenvalues for Nf_' + Nf_list[k] + ' is', sum(S[:,k])
        name = 'Nf_' + Nf_list[k]
        plt.plot(S[:,k]/sum(S[:,k]),label=str(name))

      plt.legend(loc='lower left')
      plt.yscale('log')
      plt.xlim([0,nsnap])
      plt.ylim([1e-9,1e-1]) 
      plt.savefig( path_output_figures + '/' + norm + '_sigma_filter_' + filter_name + '.png', dpi=600 )
      plt.close()
      #  plt.show()

      #for k in range(0,len(Nf_list)):
      #  name = 'Nf_' + Nf_list[k]
      #  plt.plot(S[:,k]/sum(S[:,k]),label=str(name))

      #plt.legend(loc='lower left')
      #plt.yscale('log')
      #plt.xlim([0,nsnap])
      #plt.ylim([1e-7,1]) 
      #plt.savefig('./output_figures_' + norm + '/sing_values_filter_' + filter_name + '/all_sing_values.pdf')
      #plt.close()
      ##  plt.show()

print 'End' 
  
