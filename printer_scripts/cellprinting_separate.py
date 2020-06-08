#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:09:37 2019

@author: emil
"""

'''
This script handles the printing of multiple globes. Movement of the printer is involved here.

You can print either at the same position (and change dishes in between), or print a pattern on one dish.

First adjust the tip of the needle to 25mm above the printing position in printermovement.py
'''


import serial
import time
import numpy as np
import matplotlib.pyplot as plt
import math


'''Here is what you might need to change'''



printer = serial.Serial('/dev/ttyUSB0', 115200)#you might need to change the path of the printer controller(arduino) if you re-plugged it.
pcontroller = serial.Serial('/dev/ttyUSB1', 115200,timeout=1)#you might need to change the path of the pressure controller(arduino) if you re-plugged it.

'''check which scale you are using, Satorius or GRAM, above the screen of the scale'''
scl='G' #Indicate which scale you are using #'G' or 'S'#G for GRAM, S for Sartorius

if scl=='S':#sartotius
    scale = serial.Serial('/dev/ttyACM2',9600,timeout=1, parity=serial.PARITY_EVEN)#sartotius#you might need to change the path of the printer if you re-plugged it.
elif scl=='G':
    scale = serial.Serial('/dev/ttyUSB2',baudrate=9600,timeout=None, parity=serial.PARITY_NONE)#GRAM#you might need to change the path of the printer if you re-plugged it.


density=1  #density of the bio ink in units of g/mL


'''define your printing pattern'''
'''each row of the array defines one printed globe:  [x, y, p, w]
1. the printer head moves horizontally, the movement is defined by x and y in mm.
2. the printer head moves 25mm down to the printing height, print with pressure p until mass w is printed
3. moves 25mm up back
'''
'''                [x, y, p, w]           '''
printing_pattern=[[0, 0, 3, 200],
                  [0, 0, 3, 200],
                  [0, 0, 3, 200],
                  [0, 0, 3, 200]]

'''name for output file, find them in the folder of this script'''
filename='cellprinting'










'''don't change the following unless you are sure about what you are doing'''








x_t = []
x_pm = []
x_w = []
t_start=0


plt.ioff()  #interactive mode "on" means that figures and plots remain open after program has ended
            #interactive mode "off" means that figures and plots need to be closed for the programm to continue

#sartorius scale
def getWeight_S():
    """
    Returns current weight measured by the scale, in gram. Make sure that the setting for "Print function" is set to
    "Manual" and the option "Zero Print" is set to "Off", otherwise the answer of the scale will be significantly delayed. 
    Returns:
        Floating point number of the current weight in gram. np.nan is returned in case of an error.
    """
    try:
        scale.write(b'P\r')
        answer = scale.readline()
        w=answer.decode('utf-8', 'ignore')
        weight=float(w[1:-6].replace(' ',''))
        print(weight)
        return weight
    except Exception as e:
        print ('error:',e)
        return np.nan
        #return e



#GRAM return mg
def getWeight_G():
    try:
        #scale.write(b'P\r')
        scale.reset_input_buffer()
        answer = scale.readline()
        w=answer.decode()
        weight=float(w.replace(' ','')[:-3])*1000
        print('answer=',w)
        return weight
    except Exception as e:
        print ('e:',e)
        return np.nan 





def measure():
    #call getWeight()
    if scl=='S':#sartotius
        w = getWeight_S()#measure weight
    elif scl=='G':
        w = getWeight_G()#measure weight
    #t=time.time()#debug
    global t_start
    tpstring = b' '
    #clean buffer
    pcontroller.reset_input_buffer()
    tpstring=pcontroller.readline()
    if tpstring != b' ':        
        tpstring=tpstring.decode()
        #print('tpstring = ',tpstring)#output of arduino
        t=int(tpstring.split()[0][2:])#time
        p=int(tpstring.split()[1][2:])#pressure
        if t_start == 0:
            t_start = t
        t = t - t_start
        x_t.append(t/1000)
        x_pm.append(p/(2**15)*10)
        x_w.append(w)
        #print(t)
    return 1






def move_print(x,y,p,w):
    move='G1 X{} Y{}\r'.format(x,y) 
    printer.write(move.encode('utf-8')) #move horizontally to new position
    printer.write(b'M400\r')
    printer.write(b'G1 Z-25\r') #go down 25 mm to print
    printer.write(b'M400\r')
    p=(int(p*255/10)+4).to_bytes(1, byteorder='big') #compute the 8-bit pressure string to be sent to the controller
    time.sleep(4) #wait 4 seconds until movements have finished (zero pressure)
    for i in range(5): #more waiting time, but this time measure the wight and zero pressure
        measure()
        print(i)
    if math.isnan(x_w[-1])==False:
        w_start_plat = x_w[-1] #weight before printing with new pressure values has started
    else:
        ws_id=-2
        while math.isnan(x_w[ws_id])==True:
            ws_id-=1
        w_start_plat = x_w[ws_id]
    print('w_start_plat',w_start_plat)
    pcontroller.write(p) #send target pressure to Arduino via RS232
    time.sleep(0.01) #wait 10 ms to make sure the controller has received the new pressure command and is in "send" mode
    keep_pressure = True
    while keep_pressure:
        measure()
        #print(x_w[-1]-w_start_plat)
        if x_w[-1]-w_start_plat > w: #  the actual amount is larger
            print('x_w[-1]',x_w[-1])
            keep_pressure = False
            pcontroller.write((0).to_bytes(1, byteorder='big')) #set pressure to zero
    time.sleep(0.01)
    for i in range(10):
        measure()
        print(i)
    printer.write(b'G1 Z25\r') #lift printer nozzle up by 25 mm
    printer.write(b'M400\r')
    time.sleep(10)
    
printer.write(b'G91\r') #relative coordinates
#printer.write(b'G1 Z-20\r')
printer.write(b'M400\r')
time.sleep(3)


#print following printing_pattern
for row in printing_pattern:
    move_print(row[0],row[1],row[2],row[3])


#printer.write(b'G1 Z20\r')
printer.write(b'M400\r')
move='G1 X{} Y{}\r'.format(0,0)
printer.write(move.encode('utf-8')) #move horizontally
printer.write(b'M84\r')


#handling np.nan from the scale, use the adjacent values to calculate for np.nan
error=0
for i,w in enumerate(x_w):
    if math.isnan(w)==True:
        adjacent_error=0
        k=i
        while math.isnan(x_w[k])==True:
            adjacent_error+=1
            k+=1
            error+=1
        for id_er in range(adjacent_error):
            x_w[i+id_er]=x_w[i-1]+(id_er+1)/(adjacent_error+1)*(x_w[i+adjacent_error]-x_w[i-1])
       #x_w[i]=round((x_w[i-1]+x_w[i+1])/2,1) #replace error mass value by the average of the value of the adjacent two values
       #error+=1
print('There are '+str(error)+' errors')


fig1=plt.figure(1, (10, 5))
border_width = 0.2
ax_size = [0+border_width, 0+border_width, 
           1-2*border_width, 1-2*border_width]
ax_a = fig1.add_axes(ax_size)
ax_a.plot(x_t,x_pm, '-', color = '#1f77b4', linewidth=2, label="pressure", zorder=1) #plot the  data
ax_a.set_ylabel('pressure (kPa)', color ='#1f77b4' )
ax_a.tick_params('y', colors='#1f77b4')
plt.ylim(0,np.max(x_pm*100)*1.2)
plt.legend(loc=2, numpoints=1)

ax_b = ax_a.twinx()       
ax_b.plot(x_t,x_w, '.-', color = '#d62728', linewidth=1, markerfacecolor='#d62728', markersize=6.0,markeredgewidth=0, label="mass(mg)", zorder=1) #plot the  data
ax_b.set_ylabel('mass(mg)',color = '#d62728')
ax_b.tick_params('y', colors='#d62728')    
plt.ylim(0,np.max(x_w)*1.15)
plt.legend(loc=1, numpoints=1)
plt.show()
plt.savefig(filename+'.png', dpi=300, format='png')

f = open(filename+'.txt','w')
for i in range(0,len(x_t)): 
    f.write(str(x_t[i])+'\t'+str(x_pm[i])+'\t'+str(x_w[i])+'\n')
f.close()


