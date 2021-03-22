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
__app__ = 'RadarPowerBalance'

# 
NAME_MODEL_VCO ='ZX95-2536C+'
NAME_MODEL_AMP_RF = 'ZX60-272LN+'
NAME_MODEL_ATTENUATOR ='VAT-3+'
NAME_MODEL_SPLITTER = 'ZX10-2-42'
NAME_MODEL_ANTENNA = 'Horn'
NAME_MODEL_MIXER = 'ZX05-43MH+'
NAME_MODEL_FILTER = 'MC39063'


class RadarPowerBalance(QMainWindow):
    
    def __init__(self,reference):
        
        self.radar = reference      # pointer to radarparameters
        QMainWindow.__init__(self)
        loadUi("RadarGuiAppPowers.ui",self)    # Cargs las configuraciones del diseño ! todos los atributos !
        
        self.setWindowIcon(QIcon('icon/icon_radar.ico')) 
        self.setWindowTitle("Power Balance ")
        self.setNameModel()
        self.setPowerBalance()
        #self.appRadar = selfPointer
        #self.radar = selfPointer.radar
        #self.setParameters()  
        self.ButtonExit.clicked.connect(self.RadarPowerExit)
        #self.ButtonSave.clicked.connect(self.ConfigurationSave)
        print('Init power balance')

    def RadarPowerExit(self):
        self.close()

    def setNameModel(self):
        
        self.label_vco_model.setText(NAME_MODEL_VCO)
        self.label_ampRf_model.setText(NAME_MODEL_AMP_RF)
        self.label_att_model.setText(NAME_MODEL_ATTENUATOR)
        self.label_splitter_model.setText(NAME_MODEL_SPLITTER)
        self.label_antennatx_model.setText(NAME_MODEL_ANTENNA)
        self.label_antennarx_model.setText(NAME_MODEL_ANTENNA)
        self.label_ampRf_rx_model.setText(NAME_MODEL_AMP_RF)
        self.label_mixer_model.setText(NAME_MODEL_MIXER)
        self.label_filter_model.setText(NAME_MODEL_FILTER)
        print('set models')

    def setPowerBalance(self):
        # set vco :
        self.label_vco_power.setText( str( self.radar.transmitter.vco.powerOut ) )
        # set att :
        self.label_att_gain.setText( str( -self.radar.transmitter.attenuator.gainDB ) )
        att_power_out = self.radar.transmitter.vco.powerOut - self.radar.transmitter.attenuator.gainDB
        self.label_att_power.setText( str (  att_power_out ))   
        # set amp rf tx :
        self.label__ampRf_gain.setText( '+'+str (self.radar.transmitter.ampRF.gainDB ))
        amp_tx_power_out = att_power_out + self.radar.transmitter.ampRF.gainDB 
        self.label_ampRf_power.setText( str( amp_tx_power_out) ) 
        # set splitter :
        self.label_splitter_gain.setText('+'+str(-self.radar.transmitter.splitter.gainDB ))
        splitt_power_out = amp_tx_power_out-self.radar.transmitter.splitter.gainDB
        self.label_splitter_power.setText( str( splitt_power_out ) )
        
        # set antenna tx :
        self.label_antenna_gain.setText( '+'+str(self.radar.transmitter.antenna.gainDB ) )
        self.label_antenna_gamma.setText('1')
        pire_power_out = splitt_power_out + self.radar.transmitter.antenna.gainDB
        self.label_antenna_pire.setText ( str( pire_power_out ) )

        # set medium AEL :
        ael_medium_power  = 20*np.log10(2*self.radar.medium.C/(self.radar.transmitter.vco.fmax + self.radar.transmitter.vco.fmax )) 
        ael_medium_power = ael_medium_power - 30*np.log10(4*np.pi ) - 20*np.log10( self.radar.target.distance_to_target)
        self.label_medium_ael.setText( str( np.round(ael_medium_power,2))) 

        # set target :
        self.label_target_gain.setText( str( self.radar.target.rcs.sigmadB ) )
        self.label_target_phase.setText( str('00') ) 
        self.label_target_distance.setText(str(self.radar.target.distance_to_target) ) # get it
        
        # set antenna rx :
        self.label_antenna_rx_gain.setText( str( self.radar.receiver.antenna.gainDB))
        self.label_antenna_rx_gamma.setText('1')
        power_antenna_rx =  np.round(pire_power_out  +  ael_medium_power  + self.radar.target.rcs.sigmadB + self.radar.receiver.antenna.gainDB,2 )
        self.label_antenna_rx_power.setText( str( power_antenna_rx ) )

        # set amp RF rx :
        self.label_amprRf_rx_gain.setText( '+'+str(self.radar.receiver.ampRf.gainDB)   ) 
        power_amp_rx = np.round( power_antenna_rx + self.radar.receiver.ampRf.gainDB,2)
        self.label_amprRf_rx_power.setText( str( power_amp_rx ) )   

        # set mixer :
        self.label_mixer_gain.setText( str(-self.radar.receiver.mixer.gaiDdB) )
        power_mix_rx =  power_amp_rx - self.radar.receiver.mixer.gaiDdB 
        self.label_mixer_power.setText(str(power_mix_rx) )
   
        # set filter :
        self.label_filter_gain.setText( '+'+str(self.radar.receiver.filter.gaindB ))
        power_filter = np.round( power_mix_rx + self.radar.receiver.filter.gaindB ,2)        
        self.label_filter_power.setText( str( power_filter  ) )

        # Result :
        self.label_result_power.setText( str( power_filter ))
        self.label_result_pire.setText( str( pire_power_out ))
        self.label_result_prx.setText( str(power_antenna_rx))
        #result_vpeak = np.round(np.sqrt( np.power(10,power_filter/10)*(1e-3)*2*self.radar.processor.impedance )*1e3,2 )    
        #self.label_result_vpeak.setText(str( result_vpeak ))   
        #self.label_filter_power(str())      
        #print('Set Power Balance',self.radar.transmitter.vco.powerOut )


'''
# for debug :
if __name__ == "__main__":
    AppBalancePower = QApplication(sys.argv)
    radarpower = RadarPower()
    radarpower.show()
    AppBalancePower.exec_()
'''


