#!/usr/bin/env python3
# -*- coding: utf-8 -*-




'''
This script prints at a series of pressures, with difined mass. (for rheology measurement)
'''
#choose scale
#change density for bioink
#reset printing volume and pressure
#check the name of saved file


import serial
import time
import numpy as np
import matplotlib.pyplot as plt
import math





'''Here is what you might need to change'''






pcontroller = serial.Serial('/dev/ttyUSB1', 115200,timeout=1)#you might need to change the path of the pressure controller(arduino) if you re-plugged it.

'''check which scale you are using, Satorius or GRAM, above the screen of the scale'''
scl='S' #Indicate which scale you are using #'G' or 'S'#G for GRAM, S for Sartorius
if scl=='S':#sartotius
    scale = serial.Serial('/dev/ttyACM0',9600,timeout=1, parity=serial.PARITY_EVEN)#sartotius#you might need to change the path of the printer if you re-plugged it.
elif scl=='G':
    scale = serial.Serial('/dev/ttyUSB2',baudrate=9600,timeout=None, parity=serial.PARITY_NONE)#GRAM#you might need to change the path of the printer if you re-plugged it.

density=1  #density of the bio ink in units of g/mL

'''define a increasing series of pressures (in bar) you want to use for printing'''

#p_set = np.asarray([0.15,0.3,0.45,0.6,0.75,0.9]) #in bar
p_set = np.asarray([1,1.5,2,2.5,3,3.5]) #in bar

t_min = 2 #if mass limit reached in less than t_min(second), automatically stop applying higher pressure
t_max = 40 #if mass limit not reached after t_max(second), automatically go to next pressure
w_max = 150 #mass(mg) of extruded material during one pressure plateau 

'''name for output file, find them in the folder of this script'''
filename='pressure_mass'








'''don't change the following unless you are sure about what you are doing'''







plt.close('all')
x_t = []
x_pm = []
x_pt = []
x_w = []


t_start = 0
p_target = 0



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
        x_pt.append(target_p) #only different line for cell printing
        x_pm.append(p/(2**15)*10)
        x_w.append(w)
        #print(t)
    return 1

 

pcontroller.write((0).to_bytes(1, byteorder='big') ) #send zero pressure to Arduino via RS232
target_p = 0
time.sleep(1)
measure()
while x_t[-1]<5: # wait 5 seconds
    measure()

stop_measuring = False

for i, target_p in enumerate(p_set):
    print('new target pressure is',target_p)
    p_byte = int((target_p*255/10).astype(int)+4)  #calibration
    p_byte=p_byte.to_bytes(1, byteorder='big') #convert target pressure to byte
    pcontroller.write(p_byte) #send target pressure to Arduino via RS232
    time.sleep(0.01) #wait 10 ms to make sure the controller has received the now pressure command and is in "send" mode
    scale.reset_input_buffer()
    measure()
    keep_pressure = True
    w_start_plat = x_w[-1]
    t_start_plat = x_t[-1]
    while keep_pressure:
        measure()
        if x_w[-1]-w_start_plat > w_max: #  increase since start of pressure plateau
            keep_pressure = False
            if x_t[-1]-t_start_plat < t_min: # if the threshold weight has been reached in less than 2 seconds, obviously the pressure is too high
                print('the threshold weight has been reached too quickly, obviously the pressure is too high')
                stop_measuring = True
        if x_t[-1]-t_start_plat > t_max:
            keep_pressure = False
        
    if stop_measuring:
        break
            
pcontroller.write((0).to_bytes(1, byteorder='big') ) #send zero pressure to Arduino via RS232
target_p = 0
#keep measuring for 5 more seconds
t_start_plat = x_t[-1]
while x_t[-1]-t_start_plat < 5: # 5seconds
    measure()       
 
    
'''handling np.nan from the scale, use the adjacent values to calculate for np.nan'''
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
print('There are '+str(error)+' errors, get a new scale if there are more than 5 errors')


#plot pressure&flow against time
x_w = np.array(x_w)
flow = (x_w[1:]-x_w[0:-1])/density
fig1 = plt.figure(1, (5,5))
border_w = 0.2
ax_size = [0+border_w, 0+border_w, 1-2*border_w,1-2*border_w]
ax1 = fig1.add_axes(ax_size)
plt.plot(x_t,x_pm/np.max(x_pm), '.-', label='pressure')
plt.plot(x_t[1:],flow/np.max(flow), '.-',label='weight difference')
plt.ylim(0,1.3)
#plt.xlim(10,28)
plt.legend(loc=2, numpoints=1)
plt.xlabel('time (s)')
plt.ylabel('pressure and weight difference (r.u.)')
plt.savefig(filename+'.png',dpi=300)
plt.show()

# save the measured data to a data file
f = open(filename+'.txt','w')
for i in range(0,len(x_t)): 
    f.write(str(x_t[i])+'\t'+str(x_pm[i])+'\t'+str(x_w[i])+'\t'+str(x_pt[i])+'\n')
f.close()

