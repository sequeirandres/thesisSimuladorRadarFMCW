# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
#
#  Title: Simulador de un radar de onda continua
#  Code :  SimuladorRadarFmcw.py  
#  Author: Sequeira Andres
#  gitHub: https://github.com/sequeirandres
#  Repositorio: https://github.com/sequeirandres/thesisSimuladorRadarFMCW
#
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
from PyQt5.QtWidgets import*
from PyQt5.QtWidgets import QMessageBox
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIcon 
import sys
import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

import numpy as np
import random

from scipy import signal
from scipy.signal import spectrogram

from RadarConfiguration import RadarConfiguration
from RadarPower import RadarPowerBalance

__author__ = 'sequera@andres'
__title__ = 'Radar_fmcw_simulator'

# Parametros Iniciales :

# Propiedades del amplificador de RF
AMPLIFICADOR_GAIN_DB = 10
AMPLFICADOR_BW = [1e9,2.5e9]

# Propiedades del Vco :
VCO_FMIN = 2.2e9
VCO_FMAX = 2.5e9
VCO_FCHIRP = 100
VCO_POWER_OUT_DB = 13

# Propiedades del atenuador :
ATTENUATOR_GAIN = +3
ATTENUATOR_BW = 6e9 

# Propiedades de la señal de salida :
SIGNAL_TX_PHASE = 0 

# Propiedades del splitt
SPLITT_GAIN_DB = 3 
SPLITT_NPORTS = 2 

# Propiedades de las anenas :
ANTENNA_GAIN_DB = 10 
ANTENNA_REFLEXION = 0.1 
ANTENNA_BW = [1.2e9 , 2.8e9]

# Propiedades del mixer:
MIXER_GAIN_DB = 6  # Perdidas en dB por la coversion a IF

# Parametros para el RCS (Radar Cross Section)
RCS_SIGMA = 10
RCS_PHASE = 0

# Target :
DISTANCE_TO_TARGET = 4


# Párametros del medio ambiente :
SPEED_OF_PROPAGATION = 3e8        # Velocidad de propagacion .
PERMITIVITY =  8.85418e-12   # C^{2}/N.m^{2} \epsilon_{0} e_electrico 
PERMEABILITY = 4*3.14159e-7       # mu_magnetica 
SNR = 14

# Parameters for filter analogico :
FILTER_FLOW = 10 
FILTER_FHIGH = 20e3
FILTER_GAIN_DB = 10 

# Parameters for process :
FREQ_SAMPLING = 40e3 
NFFT = 4096
SCALE = 1
PROCESS_IMPEDANCE = 50

# impedancia :
IMPEDANCIA = 50 

# Banda L :
BANDA_L_FMIN = 1e9 
BANDA_L_FMAX = 2e9 

# Banda S :
BANDA_S_FMIN = 2e9 
BANDA_S_FMAX = 4e9 

''' Define class for radar , transmitter and receiver'''

# Define class :
# Class for transmitter :
class Vco():
    def __init__(self,fmin,fmax,fchirp,powerOut):
        self.fmin = fmin
        self.fmax = fmax
        self.b = fmax-fmin 
        self.fchirp = fchirp
        self.tchirp = 1/fchirp
        self.powerOut = powerOut
        self.modulation = 'chirp'

class AmpRF():
    def __init__(self,gain,fmin,fmax):
        self.gainDB = gain
        self.fmin = fmin
        self.fmax = fmax

class Attenuator():
    def __init__(self,gainDB,BW):
        self.gainDB = gainDB
        self.BW = BW 

class Splitter():
    def __init__(self,gainDB,Nports):
        self.gainDB = gainDB 
        self.Nports = Nports

class Antenna():
    def __init__(self,gainDB,fmin, fmax):
        self.gainDB = gainDB
        self.fmin = fmin
        self.fmax = fmax

class Mixer():
    def __init__(self,gaindB):
        self.gaiDdB = gaindB
    
    def outPutIf( self , signalLo, signalRf):
        self.signal =signalLo*signalRf 

    def fftIf(self,radar):
        
        tdelay = 2*radar.target.distance_to_target/radar.medium.C

        # Retardo en frecuencia :
        fmixer = ((2*(radar.transmitter.vco.b)/radar.transmitter.vco.tchirp )*tdelay)
        # Base temporal :

        freqL = np.arange(0,16000, 0.1 )
        freqH = np.arange(radar.transmitter.vco.fmin*0.75,radar.transmitter.vco.fmin*1.25,1000000)
        signalL = np.sinc((freqL-fmixer)/20) 
        signalH = np.sinc((freqH-radar.transmitter.vco.fmin)/10000e3)

        self.freq = np.concatenate((freqL,[freqL[len(freqL)-1],freqH[0]] ,freqH , [freqH[len(freqH)-1],10e9]), axis=0)
        self.fftsignalIF = np.abs(np.concatenate((signalL, [0,0],signalH,[0,0]), axis=0))
        # add noise
        self.fftsignalIF = self.fftsignalIF +   np.random.normal( 0 , 0.005, len( self.fftsignalIF ) )

class Filter():
    def __init__(self, gaindB , flow , fhigh , fs) :
        
        fs = 40e3
        self.flow = flow 
        self.fhigh = fhigh 
        self.gaindB = gaindB
        
        b = signal.firwin(40, fhigh/fs, window=('kaiser', 8))
        w, h = signal.freqz(b)
        self.freq = (w/np.pi)*fs
        self.freq = np.concatenate(( self.freq , [fs/2,10e9] ), axis=0  )
        self.H = np.concatenate(( abs(h) , [0,0] ), axis=0 ) 


# Medium :
class Medium():
    def __init__ (self,C_luz, epsion_electic, mu_magnetic, snr  ):
        self.C = C_luz 
        self.epsilon = epsion_electic
        self.mu = mu_magnetic
        self.impedance =   np.sqrt(mu_magnetic/epsion_electic)     #obtener el campo H , en caso de necesitarlo\\ 
        self.snr = snr 
        
        powerSignal = 10e-3  # Potencia normalizada
        self.desvio =  np.sqrt(powerSignal/(np.power(10,snr/10)))                    # potencia del ruido
        self.PowerNoise = self.desvio*self.desvio
    def propagation(self,transmisor , target ):
        self.aeldB = 20*np.log10(  2*np.pi*(transmisor.vco.fmax+transmisor.vco.fmin)*target.distance_to_target*target.distance_to_target/SPEED_OF_PROPAGATION  )


# Target :

class RCS():
    def __init__(self,sigma , phase):
        self.sigma = sigma
        self.sigmadB = 10*np.log10( sigma ) 
        self.phase = phase


''' Class Target :''' 

# Target :
class Target():
    def __init__(self,distance = DISTANCE_TO_TARGET):
        self.distance_to_target = distance
        self.rcs =  RCS(RCS_SIGMA,RCS_PHASE )

''' Class Transmitter : '''

class Transmitter():
    def __init__(self):
        self.vco = Vco(VCO_FMIN,VCO_FMAX , VCO_FCHIRP , VCO_POWER_OUT_DB)
        self.attenuator = Attenuator(ATTENUATOR_GAIN,ATTENUATOR_BW)
        self.ampRF =  AmpRF(AMPLIFICADOR_GAIN_DB,BANDA_S_FMIN,BANDA_S_FMIN)
        self.splitter = Splitter(SPLITT_GAIN_DB, SPLITT_NPORTS)
        self.antenna = Antenna(ANTENNA_GAIN_DB ,BANDA_S_FMIN,BANDA_S_FMAX)

    def transmitPower(self):
        self.piredB = self.vco.powerOut-self.attenuator.gainDB-self.splitter.gainDB + self.ampRF.gainDB+self.antenna.gainDB
    
    def transmitSingnal(self):
        self.time = np.arange(0,self.vco.tchirp,1e-8)
        self.signal = np.cos(2*np.pi*( self.vco.fmin + ((self.vco.fmax-self.vco.fmin)/self.vco.tchirp)*self.time )*self.time )
        
        fs = 6000
        ts = 1/fs
        fmin = self.vco.fmin/1e6
        fmax = self.vco.fmax/1e6
        Tchrp = 10 
        fmax = (fmax-fmin)*0.5 + fmin
        NFFT = np.power(2,19)
        # Base temporal :
        time = np.arange(0,Tchrp, ts )
        # Señal chirp :
        tx_signal = np.cos(2*np.pi*( fmin + ((fmax-fmin)/Tchrp)*time )*time )
        self.freq =  np.arange(0,fs/2-fs/NFFT,fs/NFFT )*1e6
        self.fftsignal =  fft_signal = np.fft.fft(tx_signal,NFFT)
        self.fftsignal  = (np.abs(fft_signal[1:int(NFFT/2)]))/548

    def transmitPhase(self):
        self.phase = 0 


''' Class Receiver :'''

class Receiver():
    def __init__(self):
        self.antenna = Antenna(ANTENNA_GAIN_DB,BANDA_S_FMIN,BANDA_S_FMAX)
        self.ampRf =   AmpRF(AMPLIFICADOR_GAIN_DB,BANDA_S_FMIN,BANDA_S_FMIN)
        self.mixer  =  Mixer(MIXER_GAIN_DB)
        self.filter = Filter(FILTER_GAIN_DB,FILTER_FLOW,FILTER_FHIGH,FREQ_SAMPLING)

    def receivePower(self , transmisor , medium , target):
        self.powerin = transmisor.piredB - medium.aeldB + self.antenna.gainDB + target.rcs.sigmadB-10*np.log10(4*np.pi) 
        self.powerOut = self.powerin + self.antenna.gainDB + self.ampRf.gainDB - self.mixer.gaiDdB + self.filter.gaindB
    def receiveSignal(self, transmisor ,medium, target) :
        tdelay = 2*target.distance_to_target/SPEED_OF_PROPAGATION
        self.time = transmisor.time
        self.signal =  np.cos(2*np.pi*( transmisor.vco.fmin + ((transmisor.vco.b)/transmisor.vco.tchirp)*(self.time-tdelay) )*(self.time-tdelay) )
        self.signal = self.signal + np.random.normal(0, medium.desvio, len(self.time))

        ## Mejorar --
        fs = 6000
        ts = 1/fs
        fmin = transmisor.vco.fmin/1e6
        fmax = transmisor.vco.fmax/1e6
        Tchrp = 10 
        fmax = (fmax-fmin)*0.5 + fmin
        NFFT = np.power(2,19)
        # Base temporal :
        time = np.arange(0,Tchrp, ts )
        # Señal chirp :
        rx_signal = np.cos(2*np.pi*( fmin + ((fmax-fmin)/Tchrp)*time )*time ) #+ np.random.normal(0, medium.desvio, len(self.time))
        rx_signal =  rx_signal + np.random.normal(0, medium.desvio, len(time))
        self.freq =  np.arange(0,fs/2-fs/NFFT,fs/NFFT )*1e6
        self.fftsignal =  fft_signal = np.fft.fft(rx_signal,NFFT)
        self.fftsignal  = (np.abs(fft_signal[1:int(NFFT/2)]))/548

    
    def receivePhase(self, transmisor, target):
        self.phase = transmisor.phase+ target.rcs.phase

''' Class Process '''
# Class processor of radar :
class Processor() :
    def __init__(self, fchirp, fs, nfft):
        self.fchirp  = fchirp 
        self.fs = fs
        self.nfft = nfft 
        self.tchirp = 1/self.fchirp
        self.impedance = PROCESS_IMPEDANCE

    def process(self,vco,medium, target):
        self.time = np.arange(0,self.tchirp*4, 1/self.fs)
        tdelay = 2*target.distance_to_target/medium.C
        self.signal = np.sin((2*np.pi*(2*(vco.fmax-vco.fmin)/vco.tchirp )*tdelay )*self.time  ) 
        self.signal = self.signal +  np.random.normal(0, medium.desvio*10, len(self.time))
        self.freq_to_dist = medium.C*vco.tchirp/(4*vco.b)
        
        #self.scope = (self.fs/2)*medium.C*vco.tchirp/(4*vco.b)
        #self.range = (self.fs/self.nfft)*self.freq_to_dist
              # self.dist = self.freq*(  SPEED_OF_PROPAGATION *vco.tchirp/(4*(vco.fmax-vco.fmin)))
    def rangeAndScope(self,vco, medium, target) :
        self.scope = (self.fs/2)*medium.C*vco.tchirp/(4*vco.b)
        self.range = (self.fs/self.nfft)*medium.C*vco.tchirp/(4*vco.b)


    def getfft(self) :   # fft of signal and distance to target.
        self.freq =  np.arange(0,self.fs/2-self.fs/self.nfft,self.fs/self.nfft )
        self.dist = self.freq*self.freq_to_dist
        self.fftsignal =  np.fft.fft(self.signal,self.nfft)
        self.fftsignal  = np.abs(self.fftsignal[1:int(self.nfft/2)])
        self.fftsignal = self.fftsignal/np.max(self.fftsignal)   # Normalizer
        self.f,self.t,self.spectrogramSignal = spectrogram( self.signal ,self.fs , window=('tukey', 0.25), nperseg=None, noverlap=None, nfft=self.nfft, detrend='constant', return_onesided=True, scaling='density', axis=-1, mode='psd')

    def getDistanceToTarget(self) :
        print('distance to target !!') 

    def getVoltage(self,powerIn):
        self.voltage = np.sqrt(  (1e-3)*np.power(10,powerIn/10)*2*self.impedance )

    def results(self):
        positionMax = 0
        for n in range(len(self.fftsignal)):
            if  self.fftsignal[n] > 0.9  :
                positionMax = n 
                break
        self.result = Result(self.dist[positionMax],0,0,self.voltage)
    #self.rangoMin = self.dist[1]/self.scale
    #self.alcance = (self.fs/2)*SPEED_OF_PROPAGATION*vco.tchirp/(4*vco.b)


''' add a class for result '''
class Result :
    def __init__(self, distance , delta_distance , phase , voltage ) :
        self.distance  = distance 
        self.delta_distance = delta_distance
        self.phase =  phase
        self.voltage = voltage 

# Class of all :
class Radar:
    def __init__(self):
        self.transmitter = Transmitter()
        self.receiver = Receiver()
        self.medium =  Medium(SPEED_OF_PROPAGATION,PERMEABILITY,PERMEABILITY,SNR)
        self.target = Target(DISTANCE_TO_TARGET) 
        self.processor = Processor( VCO_FCHIRP ,FREQ_SAMPLING,NFFT )
        # range and scope :
        self.processor.rangeAndScope(self.transmitter.vco,self.medium ,self.target)


#radar = Radar()

# transmitir :

#radar.transmitter.transmitPower()
#radar.transmitter.transmitSingnal()
#radar.transmitter.transmitPhase()

# Propagar :
#radar.medium.propagation(radar.transmitter,radar.target)


# Recibir:
#radar.receiver.receivePower(radar.transmitter,radar.medium, radar.target )
#radar.receiver.receiveSignal(radar.transmitter, radar.medium,radar.target)
#radar.receiver.receivePhase(radar.transmitter, radar.target)

#radar.receiver.mixer.outPutIf(radar.transmitter.signal,radar.receiver.signal)

#radar.receiver.mixer.fftIf(radar)  
#radar.receiver.mixer(radar)

# Procesar:
#radar.processor.process(radar.transmitter.vco,radar.medium,radar.target)
#radar.processor.getfft()
#radar.processor.getVoltage(radar.receiver.powerOut)

'''
print('vco power Out :', radar.transmitter.vco.powerOut)
print('Atenuador :',radar.transmitter.attenuator.gainDB) 
print('Tx Amp RF Gain :', radar.transmitter.ampRF.gainDB)
print('Splitt Gain :',radar.transmitter.splitter.gainDB)
print('Tx Antena Gain :', radar.transmitter.antenna.gainDB)
print('Pire :', radar.transmitter.piredB)
print('Trasnmisor Phase :', radar.transmitter.phase)
print('ael dB', radar.medium.aeldB)
print('Power Recived :', np.round( radar.receiver.powerin ,2) )
print('power desvio :', radar.medium.desvio)
print('Scope : ', radar.processor.scope)
print('range :',radar.processor.range)
print('Power Out :',radar.receiver.powerOut)
print('Voltage @ 50 Ohms',radar.processor.voltage)

'''

''' Debug  '''

'''
plt.figure(1)
plt.plot(radar.transmitter.time[1:6000],radar.transmitter.signal[1:6000],'k') 
plt.grid(True)
plt.title('Señal Transmitida')

plt.figure(2)
plt.plot(radar.receiver.time[1:6000],radar.receiver.signal[1:6000],'k') 
plt.title('Señal Recibida')
plt.grid(True)


plt.figure(3)
plt.plot(radar.transmitter.freq,radar.transmitter.fftsignal,'k') 
plt.title('FFT signal')
plt.grid(True)


plt.figure(4)
plt.plot(radar.receiver.freq,radar.receiver.fftsignal,'k') 
plt.title('FFT signal rx')
plt.grid(True)

plt.figure(5)
plt.plot(radar.processor.time,radar.processor.signal,'k') 
plt.title('signal Io')
plt.grid(True)


plt.figure(6)
plt.plot(radar.processor.dist ,radar.processor.fftsignal,'k') 
plt.title('distance signal Io to target')
plt.grid(True)


plt.figure(10)
plt.pcolormesh(radar.processor.t,radar.processor.f*radar.processor.freq_to_dist, radar.processor.spectrogramSignal)

plt.figure(7)
plt.plot(radar.receiver.time,radar.receiver.mixer.signal,'k') 
plt.title('signal if')
plt.grid(True)



plt.figure(8)
plt.semilogx(radar.receiver.mixer.freq,radar.receiver.mixer.fftsignalIF,'k') 
plt.semilogx(radar.receiver.filter.freq,radar.receiver.filter.H, '--r')
plt.title('signal  fft If')
plt.grid(True)

plt.show()
'''

class RadarSimulador(QMainWindow):
    
    def __init__(self):
        QMainWindow.__init__(self)
        print('loading ....')
        loadUi("RadarGuiApp.ui",self)  # Load RadarGuiApp  -- configurations --
        self.initQtplot()               
        self.initMenu()                #  Load Menu - file , view , config  , about .
        self.setWindowTitle("Radar Simulator")
        self.setWindowIcon(QIcon('icon/Icono_python.ico')) 
       
        self.setRadialButton()
        self.BoxGraphics.setEnabled(False)  # Turn Off the boxedGraphics !

        self.pushButtonExit.clicked.connect(self.RadarExit)
        self.pushButtonSimular.clicked.connect(self.RadarRun)
        self.radar = Radar() # Buid a variable radar.
       
        self.setTableValues() # add scope and range
        
        # Setting scroll bar to distance : range and scope ,ready done.
        self.setScrollBar()

        self.ScrollBarTargetDistance.valueChanged['int'].connect(self.setDistanceToTarget)
        
        # Graphics options :
        self.radio_graphics_tx_time.clicked.connect( self.graphic_signal_tx_time )
        self.radio_graphics_tx_fft.clicked.connect( self.graphic_signal_tx_fft)

        self.radio_graphics_rx_time.clicked.connect( self.graphic_signal_rx_time)
        self.radio_graphics_rx_fft.clicked.connect( self.graphic_signal_rx_fft)
        
        self.radio_graphics_if_time.clicked.connect(self.graphic_signal_if_time)
        self.radio_graphics_if_fft.clicked.connect(self.graphic_signal_if_fft)

        self.radio_graphics_io_time.clicked.connect(self.graphic_signal_io_time)
        self.radio_graphics_io_fft.clicked.connect(self.graphic_signal_io_fft)
        self.radio_graphics_io_distance.clicked.connect(self.graphic_signal_io_distance)
        self.radio_graphics_io_spectrogram.clicked.connect(self.graphic_signal_io_spectrum)
        # radio_graphics_tx_time
        # radio_graphics_tx_fft 
        # self.pushButton_generate_random_signal.clicked.connect(self.update_graph)
        # self.pushButtonExit.clicked.connect(self.exit)
        # self.pushButtonSimular.clicked.connect(self.Run_Radar)
        # self.addToolBar(NavigationToolbar(self.PlotWidgetIo.canvas, self))
        # self.label_target.setText('0.0 [m]')
        # self.horizontalScrollBar.valueChanged['int'].connect(self.moveTarget)   # Muevo la barra de distancia
        # self.PlotWidgetTX.canvas.plot.subplots_adjust(left=0.05, bottom=0.150, right=0.995, top=0.9, wspace=0.2, hspace=0.2)
        # self.windowsToolConfig = windowConfigurations()
        # SET TEXTO EN LAS FIGURAS 
        
        # self.FiguraIo_ylabel.setText('Tiempo [s]')
        # self.FiguraTX_ylabel.setText('Tiempo [ms]')
        
        # self.radioButtonTime.
        # self.radioButtonIo_time.clicked.connect(self.Io_plot_time)
        # self.radioButtonIo_dist.clicked.connect(self.Io_plot_distance)
        # self.radioButtonIo_spectrograma.clicked.connect(self.Io_plot_spectrum)
        # self.radioButtonTX_time.clicked.connect(self.TX_plot_time)
        # self.radioButtonTX_fft.clicked.connect(self.TX_plot_fft)
        
      #  self.radioButtonRX_time.clicked.connect(self.RX_plot_time)
      #  self.radioButtonRX_fft.clicked.connect(self.RX_plot_fft)
        
      #  self.radioButtonIF_time.clicked.connect(self.IF_plot_time)
      #  self.radioButtonIF_fft.clicked.connect(self.IF_plot_fft)

        print('loaded ... Done ')
  # Run simulation :
    def RadarRun(self):
        print('Running simulation .....')

        # Habilitar botones de opciones graficas.
        self.BoxGraphics.setEnabled(True)  # Turn Off the boxedGraphics !            # Set Radial-button
        self.progressBar.setValue(5)
        self.setRadialButton()             # Set Radial-button
        self.progressBar.setValue(10)
        # transmitir :
        self.radar.transmitter.transmitPower()
        self.progressBar.setValue(15)
        self.radar.transmitter.transmitSingnal()
        self.progressBar.setValue(20)
        self.radar.transmitter.transmitPhase()
        # Propagar :
        self.progressBar.setValue(25)
        self.radar.medium.propagation(self.radar.transmitter,self.radar.target)
        # Recibir:
        self.progressBar.setValue(30)
        self.radar.receiver.receivePower(self.radar.transmitter,self.radar.medium, self.radar.target )
        self.progressBar.setValue(35)
        self.radar.receiver.receiveSignal(self.radar.transmitter, self.radar.medium,self.radar.target)
        self.progressBar.setValue(40)
        self.radar.receiver.receivePhase(self.radar.transmitter, self.radar.target)
        self.progressBar.setValue(45)
        self.radar.receiver.mixer.outPutIf(self.radar.transmitter.signal,self.radar.receiver.signal)
        self.progressBar.setValue(50)
        self.radar.receiver.mixer.fftIf(self.radar)  
        #radar.receiver.mixer(radar)
        self.progressBar.setValue(55)
        # Procesar:
        self.radar.processor.process(self.radar.transmitter.vco,self.radar.medium,self.radar.target)
        self.progressBar.setValue(60)
        self.radar.processor.getfft()
        self.progressBar.setValue(65)
        self.radar.processor.getVoltage(self.radar.receiver.powerOut)
        # graphics signals :
        self.progressBar.setValue(70)

        self.radar.processor.results()

        self.progressBar.setValue(80)

        self.PlotGraphicSignals() 
        # update table box :
        self.progressBar.setValue(99)
        self.updateTableValues() 

        self.progressBar.setValue(100)
        '''
        positionMax = 0
        for n in range(len(self.radar.processor.fftsignal)):
            if  self.radar.processor.fftsignal[n] > 0.9  :
                positionMax = n 
                break
        '''
        #print(self.radar.processor.result.distance )

        #self.menuBar.ViewMenu_option_balance.setEnabled(True) 
        print('simulation done !')
    
    def setScrollBar(self):
        self.ScrollBarTargetDistance.setMaximum(int(self.radar.processor.scope)*10)
        self.ScrollBarTargetDistance.setValue(self.radar.target.distance_to_target*10)
        print('set scrollbar')


    def setDistanceToTarget(self):
        self.radar.target.distance_to_target = self.ScrollBarTargetDistance.value()/10
        self.label_target_distance.setText(str( np.round(  self.radar.target.distance_to_target ,1 ) ) )
        #print(' set distance ')

#  Init radial button selection :
    def setRadialButton(self) :
        self.radio_graphics_tx_time.setChecked(True)
        self.radio_graphics_rx_time.setChecked(True)
        self.radio_graphics_if_time.setChecked(True)
        self.radio_graphics_io_distance.setChecked(True)

    def setTableValues(self):
        #print('set values')
        self.label_vco_pwr.setText(str(self.radar.transmitter.vco.powerOut))
        self.label_vco_fo.setText(str( int( (self.radar.transmitter.vco.fmax+self.radar.transmitter.vco.fmin)/2e6 )))
        self.label_vco_bw.setText(str(  int( self.radar.transmitter.vco.b/1e6 )  ))
        self.label_vco_fchirp.setText(str( self.radar.transmitter.vco.fchirp ) )
       # self.label_powers_ptx.setText(str(self.radar.vco.powerOut+self.radar.ampRf1.gainDB))
        self.label_powers_pire.setText('00.00')  # No hay potencia para tx 
        self.label_powers_ael.setText('00.00')   #
        self.label_powers_ptx.setText('00.00')
       # self.label_process_range.setText(str(self.radar.process.range ))
       # self.label_process_scope.setText(str(self.radar.process.scope))
        self.label_process_fs.setText(str( np.round( self.radar.processor.fs/1e3 ,1)))
        self.label_process_nfft.setText(str( self.radar.processor.nfft ))
        self.label_result_distance.setText('00.00 ± 00.00')
        self.label_result_phase.setText('00.00°± 00.00°')
        self.label_target_sigma.setText(str( self.radar.target.rcs.sigmadB ))
        self.label_target_phase.setText(str(self.radar.target.rcs.phase))
        self.label_target_distance.setText(str( np.round(  self.radar.target.distance_to_target ,1 ) ) )
        self.label_process_range.setText( str(  np.round(self.radar.processor.range*1e3 ,2)  ) )
        self.label_process_scope.setText( str( np.round( self.radar.processor.scope ,1))) 
        


       # self.ScrollBarTargetDistance.setMaximum(int(self.radar.processor.scope)*10)
       # self.ScrollBarTargetDistance.setValue(self.radar.target.distance_to_target*10)
        #self.radar.process.alcance = (self.radar.process.fs/2)*SPEED_OF_PROPAGATION*self.radar.transmitter.vco.tchirp/(4*self.radar.transmitter.vco.b)

    def updateTableValues(self):
        self.label_powers_pire.setText(str( np.round(self.radar.transmitter.piredB,1) ))  # No hay potencia para tx 
        self.label_result_distance.setText(str( np.round(self.radar.processor.result.distance,2) ) + '±' +  str( np.round(self.radar.processor.range,4) ))
        #self.label_result_phase.setText(str(self.radar.process.result.phase)+'°±'+'00.00°')
        self.label_powers_ael.setText(str(np.round( self.radar.medium.aeldB ,2) ))
        self.label_powers_prx.setText(str(np.round( self.radar.receiver.powerin ,2)))
        # get it from targets
        self.label_result_voltage.setText(str(self.radar.target.rcs.sigmadB ))

    #   self.label_process_range.setText( str(  np.round(self.radar.processor.range*1e3 ,2)  ) )
    #   self.label_process_scope.setText( str( np.round( self.radar.processor.scope ,1)) )

# Pannel of graphic control   ---
    def graphic_signal_tx_time(self):
        print('plot time signal tx ...  Done')
        self.PlotSignalTx('time')

    def graphic_signal_tx_fft(self):
        print('plot fft signal  tx ...  Done')
        self.PlotSignalTx('fft')

    def graphic_signal_rx_time(self):
        print('plot time signal rx ... Done')
        self.PlotSignalRx('time')

    def graphic_signal_rx_fft(self):
        print('plot fft  signal fft ... Done')
        self.PlotSignalRx('fft')

    def graphic_signal_if_time(self):
        print('plot time signal IF  ... Done')
        self.PlotSignalIf('time')

    def graphic_signal_if_fft(self):
        print('plot fft signal IF  ... Done')
        self.PlotSignalIf('fft')

    def graphic_signal_io_time(self):
        print('plot time signal Io ... Done')
        self.PlotSignalIo('time')
    
    def graphic_signal_io_fft(self):
        print('plot fft signal Io ... Done')
        self.PlotSignalIo('fft')
    
    def graphic_signal_io_distance(self):
        print('plot distance signal Io ...  Done')
        self.PlotSignalIo('distance')

    def graphic_signal_io_spectrum(self):
        print('plot spectrum signal Io ... Done')
        self.PlotSignalIo('spectrum')

# ---- finalizar el programa ---------

    def RadarExit(self):
        print('End Of App Radar ')
        sys.exit()

# ---- Inicializar los graficos -------

    def initQtplot(self):

        self.GraphicSignalTx.canvas.axes.clear()
        self.GraphicSignalTx.canvas.axes.grid(True)    
        self.GraphicSignalTx.canvas.draw()

        self.GraphicSignalRx.canvas.axes.clear()
        self.GraphicSignalRx.canvas.axes.grid(True) 
        self.GraphicSignalRx.canvas.draw()

        self.GraphicSignalIf.canvas.axes.clear()
        self.GraphicSignalIf.canvas.axes.grid(True) 
        self.GraphicSignalIf.canvas.draw()

        self.GraphicSignalIo.canvas.axes.clear()
        self.GraphicSignalIo.canvas.axes.grid(True) 
        self.GraphicSignalIo.canvas.draw()

#  ---- graficos ------

    def PlotGraphicSignals(self):
        print('Init Graphics ')
        self.PlotSignalTx('time')
        self.PlotSignalRx('time')
        self.PlotSignalIf('time')
        self.PlotSignalIo('distance')
        print('Graphics signals .... Done')

    def PlotSignalTx(self,option='time'):
        self.GraphicSignalTx.canvas.axes.clear()

        if option == 'time' :
            self.GraphicSignalTx.canvas.axes.plot(self.radar.transmitter.time[1:6000]*1e6, self.radar.transmitter.signal[1:6000],'k')
            self.lable_graphics_tx.setText('Time [ns]')
            self.GraphicSignalTx.canvas.axes.set_ylabel('Amplitude [v]')
            #self.PlotWidgetTX.canvas.axes.set_xlabel('Frecuencia [Hz]')   # AGREGO XLABEL 
            #self.PlotWidgetTX.canvas.axes.set_ylabel('Amplitud [V]')   # AGREGO YLABEL 
        elif option == 'fft':
            self.GraphicSignalTx.canvas.axes.plot(self.radar.transmitter.freq , self.radar.transmitter.fftsignal,'k')
            self.GraphicSignalTx.canvas.axes.set_xlim([ self.radar.transmitter.vco.fmin-300e6 ,  self.radar.transmitter.vco.fmax+300e6 ])
            self.lable_graphics_tx.setText('Frequency [GHz]')
            self.GraphicSignalTx.canvas.axes.set_ylabel('| Amplitude |')
            #self.GraphicSignalTx.canvas.axes.set_xscale('log')
        else:
            pass
        
        self.GraphicSignalTx.canvas.axes.grid(True)
        self.GraphicSignalTx.canvas.draw()
        
        #self.PlotWidgetTX.subplots_adjust(left=0.05, bottom=0.150, right=0.995, top=0.9, wspace=0.2, hspace=0.2)
        #self.plotTxGrafico.canvas.axes.plot.subplots_adjust(left=0.05, bottom=0.150, right=0.995, top=0.9, wspace=0.2, hspace=0.2)
        # self.plotTxGrafico.tight_layout() # tight_layout()       #     plt.tight_layout()

    def PlotSignalRx(self,option='time'):
        self.GraphicSignalRx.canvas.axes.clear()

        if option == 'time' :
            self.GraphicSignalRx.canvas.axes.plot(self.radar.receiver.time[1:6000]*1e6, self.radar.receiver.signal[1:6000] ,'k')
            self.lable_graphics_rx.setText('Time [ns]')
            self.GraphicSignalRx.canvas.axes.set_ylabel('Amplitude [v]')
            #self.PlotWidgetTX.canvas.axes.set_xlabel('Frecuencia [Hz]')   # AGREGO XLABEL 
            #self.PlotWidgetTX.canvas.axes.set_ylabel('Amplitud [V]')   # AGREGO YLABEL 
        elif option == 'fft':
            self.GraphicSignalRx.canvas.axes.plot(self.radar.receiver.freq , self.radar.receiver.fftsignal,'k')
            self.GraphicSignalRx.canvas.axes.set_xlim([ self.radar.transmitter.vco.fmin-300e6 ,  self.radar.transmitter.vco.fmax+300e6 ])
            self.lable_graphics_rx.setText('Frequency [GHz]')
            self.GraphicSignalRx.canvas.axes.set_ylabel('| Amplitude |')
            #self.GraphicSignalRx.canvas.axes.set_xscale('log')
        else:
            pass

        self.GraphicSignalRx.canvas.axes.grid(True)
        self.GraphicSignalRx.canvas.draw()

    def PlotSignalIf(self,option='time'):
        self.GraphicSignalIf.canvas.axes.clear()
        #self.PlotWidgetIF.canvas.axes.plot(self.t,self.signal_3 ,'black')
        if option == 'time' :
            self.GraphicSignalIf.canvas.axes.plot(self.radar.receiver.time*1e3 , self.radar.receiver.mixer.signal,'k')
            self.GraphicSignalIf.canvas.axes.set_xlim([0,self.radar.transmitter.vco.tchirp*1e3/2])
            self.lable_graphics_if.setText('Time [ms]') 
            self.GraphicSignalIf.canvas.axes.set_ylabel('Amplitude [v]')
                       #self.PlotWidgetTX.canvas.axes.set_xlabel('Frecuencia [Hz]')   # AGREGO XLABEL 
            #self.PlotWidgetTX.canvas.axes.set_ylabel('Amplitud [V]')   # AGREGO YLABEL 
        elif option == 'fft':
            self.GraphicSignalIf.canvas.axes.plot(self.radar.receiver.mixer.freq, self.radar.receiver.mixer.fftsignalIF,'k')
            self.GraphicSignalIf.canvas.axes.plot(self.radar.receiver.filter.freq, self.radar.receiver.filter.H ,'--r')
            self.GraphicSignalIf.canvas.axes.set_xscale('log')
            self.lable_graphics_if.setText('Frequency [Hz]') 
            self.GraphicSignalIf.canvas.axes.set_ylabel('| Amplitude |')
        else:
            pass

        self.GraphicSignalIf.canvas.axes.grid(True)
        self.GraphicSignalIf.canvas.draw()

    def PlotSignalIo(self,option='distance'):
        self.GraphicSignalIo.canvas.axes.clear()
        if option == 'time' :
            self.GraphicSignalIo.canvas.axes.plot(self.radar.processor.time,self.radar.processor.signal ,'k')
            self.GraphicSignalIo.canvas.axes.set_xlim([0,self.radar.transmitter.vco.tchirp/2])
            self.lable_graphics_io.setText('Time [s]')
            self.GraphicSignalIo.canvas.axes.set_ylabel('Amplitude [v]')
            #self.PlotWidgetTX.canvas.axes.set_xlabel('Frecuencia [Hz]')   # AGREGO XLABEL 
            #self.PlotWidgetTX.canvas.axes.set_ylabel('Amplitud [V]')   # AGREGO YLABEL 
        elif option == 'fft':
            self.GraphicSignalIo.canvas.axes.plot(self.radar.processor.freq ,self.radar.processor.fftsignal,'k')
            if self.radar.processor.result.distance < self.radar.processor.scope/2 :
                self.GraphicSignalIo.canvas.axes.set_xlim(0,self.radar.processor.fs/4)
            else :
                self.GraphicSignalIo.canvas.axes.set_xlim(self.radar.processor.fs/4,self.radar.processor.fs/2)
            self.lable_graphics_io.setText('Frequency [Hz]')
            self.GraphicSignalIo.canvas.axes.set_ylabel('| Amplitud |')
        elif option == 'distance':
            self.GraphicSignalIo.canvas.axes.plot(self.radar.processor.dist ,self.radar.processor.fftsignal,'k')
            if self.radar.processor.result.distance < self.radar.processor.scope/2 :
                self.GraphicSignalIo.canvas.axes.set_xlim(0,self.radar.processor.scope/2)
            else :
                self.GraphicSignalIo.canvas.axes.set_xlim(self.radar.processor.scope/2,self.radar.processor.scope)
            self.lable_graphics_io.setText('Distance [m]') 
            self.GraphicSignalIo.canvas.axes.set_ylabel('| Amplitude |')
            #self.FiguraIo_ylabel.setText('Distanacia [m]')
        elif option == 'spectrum':   #colormesh
            self.GraphicSignalIo.canvas.axes.contourf(self.radar.processor.t ,self.radar.processor.f*self.radar.processor.freq_to_dist,self.radar.processor.spectrogramSignal,cmap='jet')
            self.GraphicSignalIo.canvas.axes.set_ylim(self.radar.processor.result.distance-1.5,self.radar.processor.result.distance+1.5)
            self.GraphicSignalIo.canvas.axes.set_ylabel('Distance [m]')   #xlabel('Distance[m]')
            #self.GraphicSignalIo.canvas.axes.colorbar()
            self.lable_graphics_io.setText('Time [s]') 
        else :
            pass

        self.GraphicSignalIo.canvas.axes.grid(True)
        self.GraphicSignalIo.canvas.draw()

# -------  BARRA DE MENU DEL PROGRAMA -----------

    def initMenu(self):

        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        editMenu = mainMenu.addMenu('Edit')
        viewMenu = mainMenu.addMenu('View')
        searchMenu = mainMenu.addMenu('Search')
        toolsMenu = mainMenu.addMenu('Tools')
        helpMenu = mainMenu.addMenu('Help')
        

        # Atributos para la Nueva barra de menu
        fileMenu_Option_Exit = QAction('Exit', self)
        fileMenu_Option_New = QAction('New',self)
        fileMenu_Option_Open = QAction('Open',self)
        fileMenu_Option_Save = QAction('Save',self)

        fileMenu_Option_Exit.triggered.connect(self.RadarExit)    # Cerrar la aplicacion
    
        helpMenu_Option_Info = QAction('About',self)
        helpMenu_Option_Info.triggered.connect(self.About)        
        
        toolsMenu_option_config = QAction('setting',self)
        toolsMenu_option_config.triggered.connect(self.toolsMenuConfigurations)

        ViewMenu_option_balance = QAction('power balance',self)
        #ViewMenu_option_balance.setEnabled(False)  # No enable until you simulated

        ViewMenu_option_balance.triggered.connect(self.viewMenuBalancePower)

        # Submenu para el 1er menu fileMenu 
        fileMenu.addAction(fileMenu_Option_New)
        fileMenu.addAction(fileMenu_Option_Open)
        fileMenu.addAction(fileMenu_Option_Save)
        fileMenu.addAction(fileMenu_Option_Exit)
        
        toolsMenu.addAction(toolsMenu_option_config)

        # Submenu para el 4to menu HelpMenu
        helpMenu.addAction(helpMenu_Option_Info)

        
        # Submenu for View Power Balance :
        viewMenu.addAction(ViewMenu_option_balance)


        #fileMenu.addMenu()
        #fileMenu.addMenu("Abrir")
        #fileMenu.addMenu("Guardar")
        #fileMenu.addMenu("Exit")
        #exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
       # exitButton.setShortcut('Ctrl+Q')
        #exitButton.setStatusTip('Exit application')
        #exitButton.triggered.connect(self.exit)    # Cerrar la aplicacion

    def toolsMenuConfigurations(self):
        self.windowsToolConfig = RadarConfiguration(self)  # send self.radar
        self.windowsToolConfig.show()


    def viewMenuBalancePower(self) :
        self.windowAppPowerBalance = RadarPowerBalance(self.radar) 
        self.windowAppPowerBalance.show()


    def About(self) :
        print('About')
        TextAbout = 'Radar Created by @andres' + '\n' + 'Simulador ver-01.2020'
        QMessageBox.question(self, 'Radar-About', TextAbout, QMessageBox.Ok, QMessageBox.Ok)

#  class of configurations :
'''
class windowConfigurations(QMainWindow):
     
    def __init__(self):
        
        QMainWindow.__init__(self)
        loadUi("RadarGuiAppConfiguration.ui",self)    # Cargs las configuraciones del diseño ! todos los atributos !
        self.pushButtonConfigExit.clicked.connect(self.ConfigExit)
        print('Init Configuration  ... Done !')
        
    def ConfigExit(self):
        print('Ennd Configuration  ...  Done !')
        self.close()
'''
#  __main__ :

if __name__ == "__main__":
    app = QApplication(sys.argv)
    Radar = RadarSimulador()
    Radar.show()
    app.exec_()

