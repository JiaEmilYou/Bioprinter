#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:12:19 2019

@author: emil
"""

'''
This script controls printing with defined pressure and time
'''

import serial
import time
printer = serial.Serial('/dev/ttyUSB0', 115200) #you might need to change the path of the printer if you re-plugged it.
pcontroller = serial.Serial('/dev/ttyUSB1', 115200,timeout=1) #you might need to change the path of the pressurecontroller(arduino) if you re-plugged it.



'''This is you input'''
p=4 #define pressure in bar, the measure pressure can be found in console during printing
t=2 #define time in second (at least 0.02 s)


'''don't change the following unless you are sure about what you are doing'''
p=(int(p*255/10)+4).to_bytes(1, byteorder='big')
pcontroller.write(p)
time.sleep(t-0.01)

pcontroller.reset_input_buffer()
tpstring=pcontroller.readline()
if tpstring != b' ':        
    tpstring=tpstring.decode()
    #print('tpstring = ',tpstring)#output of arduino
    t=int(tpstring.split()[0][2:])
    p=int(tpstring.split()[1][2:])
    print('measured pressure:',p/(2**15)*10)
time.sleep(0.01)
p=0
p=(int(p*255/10)+4).to_bytes(1, byteorder='big')
pcontroller.write(p)
