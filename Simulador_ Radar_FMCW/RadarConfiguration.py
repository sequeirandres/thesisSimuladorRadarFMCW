from PyQt5.QtWidgets import*
from PyQt5.QtWidgets import QMessageBox
from PyQt5.uic import loadUi
from PyQt5.QtGui import QIcon 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

import numpy as np
import random
import sys

from scipy import signal

__author__ = 'andres@sequeira'
__app__ = 'RadarConfugrations'

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
MIXER_GAIN_DB = 3  # Perdidas en dB por la coversion a IF

# Parametros para el RCS (Radar Cross Section)
RCS_SIGMA = 10
RCS_PHASE = 0

# Target :
DISTANCE_TO_TARGET = 4


# Párametros del medio ambiente :
SPEED_OF_PROPAGATION = 3e8        # Velocidad de propagacion .
PERMITIVITY =  8.85418e-12   # C^{2}/N.m^{2} \epsilon_{0} e_electrico 
PERMEABILITY = 4*3.14159e-7       # mu_magnetica 
SNR = 10

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
# -----------------------------------------------

class RadarConfiguration(QMainWindow):
    
    def __init__(self,selfPointer):
        
        QMainWindow.__init__(self)
        
        loadUi("RadarGuiAppConfiguration.ui",self)    # Cargs las configuraciones del diseño ! todos los atributos !
        self.setWindowIcon(QIcon('icon/Icono_setting.ico')) 
        self.setWindowTitle("Setting ")
        self.appRadar = selfPointer
        self.radar = selfPointer.radar
        self.setParameters()  
        self.ButtonConfigExit.clicked.connect(self.RadarConfigurationExit)
        self.ButtonSave.clicked.connect(self.ConfigurationSave)
        print('Init Configuration  ... Done !')
        

    '''  
    def __init__(self,sarValues):
        
        QMainWindow.__init__(self)
        loadUi("RadarGuiAppConfiguration.ui",self)    # Cargas las configuraciones del diseño ! todos los atributos !
        #self.loadParameters(sarValues)
        self.sar = sarValues  
        self.setParameters()               
        self.setWindowTitle("Configuraciones ")
        self.setWindowIcon(QIcon('Icono_config.ico')) 
        self.pushButton_Exit.clicked.connect(self.ConfigurationExit)
        self.pushButton_Save.clicked.connect(self.ConfigurationSave)
        self.pushButton_Default.clicked.connect(self.ConfigurationDefault)  
        self.pushButton_value.clicked.connect(self.printValues)                      
        print('loading configuration .... done!')
    '''
#    def loadParameters(self,sarValues) :
#        self.f_max = sarValues.f_max
#        self.f_min = sarValues.f_min
#        self.f_chrp = sarValues.f_chrp
#        self.t_chrp = sarValues.t_chrp
#        self.adc_fsample = sarValues.adc_fsample
#        self.adc_nfft = sarValues.adc_nfft
#        self.tx_potencia_dbm = sarValues.tx_potencia_dbm
#        self.tx_srn = sarValues.tx_srn

    def setParameters(self):

        self.BoxVco_fmin.setValue(self.radar.transmitter.vco.fmin/1e9)
        self.BoxVco_fmax.setValue(self.radar.transmitter.vco.fmax/1e9) 
        self.BoxVco_fchirp.setValue(self.radar.transmitter.vco.fchirp )
        self.BoxVco_power.setValue(self.radar.transmitter.vco.powerOut ) 
        self.BoxVco_phase.setValue(0 )

        self.BoxAttenuator_gain.setValue( self.radar.transmitter.attenuator.gainDB )

        self.BoxSplitter_gain.setValue( self.radar.transmitter.splitter.gainDB )
        self.BoxSplitter_nport.setValue( self.radar.transmitter.splitter.Nports)

        self.BoxAmpRF_gain.setValue( self.radar.transmitter.ampRF.gainDB)
        self.BoxAmpRF_fmin.setValue( self.radar.transmitter.ampRF.fmin )
        self.BoxAmpRF_fmax.setValue( self.radar.transmitter.ampRF.fmax )

        self.BoxMixer_gain.setValue ( self.radar.receiver.mixer.gaiDdB )

        self.BoxProcess_fs.setValue ( self.radar.processor.fs/1e3 )
        self.BoxProcess_nfft.setValue( self.radar.processor.nfft )
      
        self.BoxAntenna_gain.setValue(self.radar.transmitter.antenna.gainDB)  
      
        self.BoxFilter_gain.setValue( self.radar.receiver.filter.gaindB)
        self.BoxFilter_flow.setValue( 0)
        self.BoxFilter_fhigh.setValue( 15)
      
        self.BoxTarget_sigma.setValue(self.radar.target.rcs.sigmadB)
        self.BoxTarget_phase.setValue( self.radar.target.rcs.phase)
        self.BoxTarget_distance.setValue( self.radar.target.distance_to_target)


       # self.doubleSpinBox_f_min.setValue(self.sar.f_min/1e9)
       # self.doubleSpinBox_fs.setValue(self.sar.adc_fsample/1e3)
       # self.spinBox_nfft.setValue(self.sar.adc_nfft)
       # self.doubleSpinBox_f_chrp.setValue(self.sar.f_chrp)
       # self.doubleSpinBox_TXpotencia.setValue(self.sar.tx_potencia_dbm)
       # self.doubleSpinBox_SNR.setValue(self.sar.tx_srn)
       # self.spinBox_scale.setValue(self.sar.scale)

    def validateValues(self):
        
        if(self.BoxVco_fmax.value() <= self.BoxVco_fmin.value()):
            QMessageBox.about(self, "Error ", "F max <= F min ")
            return False
        else :
            return True

    def ConfigurationSave(self):

        if(self.validateValues()):
            print('Values Valid ')
            #self.pushButton_Default.setEnabled(False)
            #self.groupSetting.setEnabled(False)
            # Upgrade Values
            # vco :
            self.radar.transmitter.vco.fmin = int (np.round(self.BoxVco_fmin.value()*1e9,0))
            self.radar.transmitter.vco.fmax = int(np.round (self.BoxVco_fmax.value()*1e9 ,0))
            self.radar.transmitter.vco.b = self.radar.transmitter.vco.fmax-self.radar.transmitter.vco.fmin
            self.radar.transmitter.vco.fchirp = int( np.round( self.BoxVco_fchirp.value() ,0)) 
            self.radar.transmitter.vco.tchirp = 1/self.radar.transmitter.vco.fchirp
            self.radar.transmitter.vco.powerOut = self.BoxVco_power.value()

            # attenuator: 
            self.radar.transmitter.attenuator.gainDB =  self.BoxAttenuator_gain.value()
            
            # Amp RF: 
            self.radar.transmitter.ampRF.gainDB = self.BoxAmpRF_gain.value()
            self.radar.transmitter.ampRF.fmin = self.BoxAmpRF_fmin.value()
            self.radar.transmitter.ampRF.fmax = self.BoxAmpRF_fmax.value()

            # Process :
            self.radar.processor.fs =  int(self.BoxProcess_fs.value()*1e3)
            self.radar.processor.nfft = int(self.BoxProcess_nfft.value())

            #self.sar.scope = (self.sar.adc_fsample/2)*(CLUZ*self.sar.t_chrp)/(4*self.sar.B)
            #self.sar.scale = self.spinBox_scale.value()
            #self.sar.scope = self.sar.scope*self.sar.scale

            # Upgrade Labels
            #self.sar.label_f_max_value.setText( str( self.sar.f_max/1e9 ))
            #self.sar.label_f_min_value.setText( str( self.sar.f_min /1e9))
            #self.sar.label_f_chirp_value.setText( str( self.sar.f_chrp ))    
            #self.sar.label_fs_value.setText( str(self.sar.adc_fsample/1e3 )) 
            #self.sar.label_nfft_value.setText(str(self.sar.adc_nfft))
            #self.sar.label_bw_value.setText(str(self.sar.B/1e6))
            #self.sar.label_scope_value.setText(str(self.sar.scope))
            #self.sar.label_scale_value.setText(str(self.sar.scale))


            print('Save Configuration  ...  Done! ')
            self.printValues()
            # Actualizo el tablero :
            self.appRadar.setTableValues()

        else :
            print('Error - Values No valid ')


    def ConfigurationDefault(self):
        
        self.loadParametersDefault()
        self.setParameters()      
        print('Default ... Done!')

    def RadarConfigurationExit(self):
        
        print('Exit Configuration ...  Done! ')
        self.close()


    def printValues(self):
        
        print('-----------------------------------------')
        print('        - Parametros del Radar -         ')
        print('-----------------------------------------')
        print(' Frecuencia maxima   [GHz] : ',self.radar.transmitter.vco.fmax) 
        print(' Frecuencia Minima   [GHz] : ',self.radar.transmitter.vco.fmin)
        print(' Ancho de Banda      [MHz] : ',self.radar.transmitter.vco.b)
        print(' Frecuencia de Chirp  [Hz] : ',self.radar.transmitter.vco.fchirp)
        print(' Periodo de chirp      [S] : ',self.radar.transmitter.vco.tchirp)
        print(' Frecuencia de sample [KHz]: ',self.radar.processor.fs)
        print(' Puntos para la FFT        : ' ,self.radar.processor.nfft)
       # print(' Alcance Max               : ',self.radar.processor.scope)
       # print(' Rango  min                : ',self.radar.processor.range)
       

    def loadParametersDefault(self):
        self.f_max = F_MIN
        self.f_min = F_MIN
        self.f_chrp = F_CHRP
        self.t_chrp = T_CHRP  #sarValues.t_chrp
        self.adc_fsample =  FS  #sarValues.adc_fsample
        self.adc_nfft =  NFFT #sarValues.adc_nfft
        self.tx_potencia_dbm =  TX_DBM #sarValues.tx_potencia_dbm
        self.tx_srn =  SNR  #sarValues.tx_srn

    