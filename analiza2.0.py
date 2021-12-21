#!/usr/bin/python
from scipy import signal
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.animation as animation
from scipy.stats import norm
from matplotlib.gridspec import GridSpec
from random import randint







def update_line(num, data, line):
    i = randint(0,len(data[0])-1024)
    #line.set_data(data[..., :num])
    line.set_data(data[0,0:1024],data[1,i:i+1024])
     # En este caso equivalente a:  line.set_data(data[:, :num])
    return line,

def update_linepsd(num,data, line):
    i = randint(1,10)
    #line.set_data(data[..., :num])
    line.set_data(data[0,0:512],data[1,i*512+i:(i+1)*512+i])
     # En este caso equivalente a:  line.set_data(data[:, :num])
    return line,

########################### Archivos ##########################################
#filetxt =  "Medicion2_op_normal_secundario_RA6_19agosto_2015.txt"   
filetxt = "Medicion_op_normal_secundario_RA6.txt"
#filetxt = "Medicion_op_normal_secundario_RA6_13_diciembre_2017.txt"
###############################################################################


########################### Datos iniciales para procesamiento ################
fs = 100
nfft = 1024
nbin=40
j=0
###############################################################################

########################### Lectura de archivo ################################
x = pd.read_csv(filetxt, sep="\t", names=["tiempo","voltage"]) #
############################################################################### 

######################### Histograma ##########################################
x_mean = x.voltage - np.mean(x.voltage)
mu = 0
sigma = np.std(x_mean)
kurt = round(sp.stats.kurtosis(x_mean),2)
skew = round(sp.stats.skew(x_mean),2)
###############################################################################    

######################## Textos ###############################################
tex0 = "Analizador de señales\n"+filetxt.replace(".txt","").replace("_"," ")
tex1 = "- Kurtosis = " + str(kurt) +"\n"+ "- Skewness = "+ str(skew)   
tex2 = "PSD con Fs= "+str(fs)+" Hz"
tex3 = "Método de Welch "+"\n"+"Ndatos="+str(len(x.tiempo))+"\n"+"Nfft ="+str(nfft)   
###############################################################################

#x0 = [x.tiempo[j*1024:(j+1)*1024],x.voltage[j*1024:(j+1)*1024]] 
#x.tiempo[j*1024:(j+1)*1024], x.voltage[j*1024:(j+1)*1024]]
#data = np.c_[x.tiempo[j*1024:(j+1)*1024],x.voltage[j*1024:(j+1)*1024]].T

######################### Ploteo de la señal cruda ############################
data = np.c_[x.tiempo, x.voltage].T
gs = GridSpec(3, 2)
fig1 = plt.figure(1)   
plt.subplot(gs[0,:])         
l, = plt.plot([],[],'g', label = 'Señal adquirida', animated='True')

line_ani  = animation.FuncAnimation(fig1, update_line, 1024, fargs=(data, l),
                                   interval=500, blit=True)
#plt.plot(x.tiempo, x.voltage, label = 'Señal adquirida')
plt.legend(loc="upper right")
plt.title(tex0)
plt.xlim(0, 10.24)
plt.ylim(4, 5)
plt.xlabel('Tiempo [s]')
plt.ylabel('Voltaje [v]')
plt.grid() 
###############################################################################    

############################### Ploteo de estadística #########################
plt.subplot(gs[1,:-1])
#plt.figure()
n, bins, patches = plt.hist(x_mean,nbin,normed=1, facecolor='green', alpha=0.5,label = "Histograma")
gaus_p = norm.pdf(bins,mu,sigma) # densidad de probabilidad de una distribución normal con ese mu y sigma
plt.plot(bins, gaus_p, 'r--', label = "Ajuste Normal")
plt.legend(loc="upper right", title=tex1)
#plt.title("Histograma de las fluctuaciones")
plt.xlabel("Ruido [V]")
plt.ylabel('Probability')    
## densidad acumulada
His, X = np.histogram(x_mean, bins=nbin, normed = True, density = True) #En His está el valor de frecuencia normado y en X está el ancho de bins
dx = X[1] - X[0] #
y_acum = np.cumsum(His)*dx #Alto x ancho es el valor de probabilidad
gaus_a = norm.cdf(bins, mu, sigma)
#tex_acum = "Prob acumuluda"
plt.subplot(gs[2,:-1])
plt.plot(X[1:], y_acum, color='green', alpha=0.5, label = "Prob acumulada")
plt.xlabel("Ruido [V]")
plt.ylabel('Cumuled Probability')
plt.plot(bins, gaus_a, 'r--', label = "Ajuste Normal")
plt.legend(loc="lower right")
###############################################################################
 
####################### cálculo de la densidad espectral ######################
f, Pxx_den = signal.welch(x.voltage, fs, nperseg=nfft)
for j in range(1,11):
    f0, p0 = signal.welch(x.voltage[(j-1)*10000:j*10000], fs, nperseg=nfft)
    p0 = signal.savgol_filter(p0,5,3) #filtrado smooth
    #Pxx_den = np.c_[Pxx_den, p0]
    Pxx_den = np.concatenate([Pxx_den, p0])
    #f = np.c_[f,f0]
    f = np.concatenate([f,f0])

datapsd = np.c_[f,Pxx_den].T    
plt.subplot(gs[1:3,1])
#fig2 = plt.figure(2)
l2, = plt.loglog([],[],'g', label = 'Power Spectral Density', animated='True')

line_ani2  = animation.FuncAnimation(fig1, update_linepsd, 512, fargs=(datapsd, l2),
                                   interval=500, blit=True)
#plt.loglog(f, Pxx_den,'g',label = tex2)
#plt.loglog(f[:20], Pxx_den[:20],'go')
plt.xlim(0.08,50)
plt.ylim(4E-6,0.005)
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [V**2/Hz]')
plt.grid()
plt.legend(loc="upper right", title= tex3)
plt.show()
###############################################################################