#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 10:39:09 2019

@author: emil
"""

'''
This script controls the movement of the printer.
'''


import serial
printer = serial.Serial('/dev/ttyUSB0', 115200, timeout = 1)#you might need to change the path of the printer if you re-plugged it.




'''This is you input'''
movement=[0,0,29]#Input movement as x,y,z coordinates, in mm.






'''don't change the following unless you are sure about what you are doing'''






x,y,z=movement
move_string="G0 X{0} Y{1} Z{2}  E300 F1000\r".format(x+100,y+100,z+100)#the origin point is 'X100 Y100 Z100'
move_byte=str.encode(move_string)
inputstring = ''
while ("wait" not in str(inputstring)) and ("fail" not in str(inputstring)):#This loop continuous as long as the printer is busy
    while printer.inWaiting():#ask whether there is a input from the printer, if not, go to next loop to ask pcontroller
        inputstring = printer.read(printer.inWaiting()) #'wait' in string(inputstring)
        print(inputstring)
printer.write(b'G92 X100 Y100 Z100\r') #origin point
printer.write(b'G90\r')
printer.write(b'M203 Z300 \r')
printer.write(move_byte) 
printer.write(b'M84\r')
