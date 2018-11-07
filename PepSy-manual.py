'''
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PepSy - An open-source peptide synthesizer
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Developed by:
Dr. Hariprasad Gali, Ph.D.
Associate Professor of Research
Department of Pharmaceutical Sciences
College of Pharmacy
The University of Oklahoma Health Sciences Center
Oklahoma City, OK 73117

Email: hgali@ouhsc.edu.

Tested only with Python 3.5.0
First full version - December 16, 2015
Last update - November 06, 2018

This script is written for the manual control of PepSy and cleaning amino acid/reagent line.
Save device configuration file (config.txt) in the same folder where PepSy.py and PepSy-manual.py scripts are saved.
Only Arduino UNO digital pins are used.
COM port numbers of VICI stream selector valve (ps) and Arduino UNO (board) needs to be updated in the device configuration file according to their current assignment on the PC.

VICI CHEMINERT low pressure 24 position stream selector valve
Position 1 - Air
Position 2 - DMF
Position 3 - DCM
Position 4 - Piperidine solution
Position 5 - DIPEA solution
Position 6 - HOBT solution
Position 7 - HBTU solution
Position 8 to 24 - Amino acid solutions or other reagent solutions
'''

import time
import serial
import _thread
import configparser
from tkinter import *
from pyfirmata import Arduino, util

# Functions
def n2On():
    n2OnBtn.config(state=DISABLED)
    n2OffBtn.config(state=NORMAL)
    n2.write(1)

def n2Off():
    n2OffBtn.config(state=DISABLED)
    n2OnBtn.config(state=NORMAL)
    n2.write(0)

def ventOn():
    ventOnBtn.config(state=DISABLED)
    ventOffBtn.config(state=NORMAL)
    vent.write(1)

def ventOff():
    ventOffBtn.config(state=DISABLED)
    ventOnBtn.config(state=NORMAL)
    vent.write(0)

def reagentOn():
    reagentOnBtn.config(state=DISABLED)
    reagentOffBtn.config(state=NORMAL)
    reagent.write(1)

def reagentOff():
    reagentOffBtn.config(state=DISABLED)
    reagentOnBtn.config(state=NORMAL)
    reagent.write(0)

def wasteOn():
    wasteOnBtn.config(state=DISABLED)
    wasteOffBtn.config(state=NORMAL)
    waste.write(1)

def wasteOff():
    wasteOffBtn.config(state=DISABLED)
    wasteOnBtn.config(state=NORMAL)
    waste.write(0)

def primeOn():
    primeOnBtn.config(state=DISABLED)
    primeOffBtn.config(state=NORMAL)
    prime.write(1)

def primeOff():
    primeOffBtn.config(state=DISABLED)
    primeOnBtn.config(state=NORMAL)
    prime.write(0)

def pumpon(v): # v is volume (integer) to be pumped in microliters, pumping rate is 1.2 ml/min, this function is required as a solenoid valve based micro pump is being used
    pumpOnBtn.config(state=DISABLED)
    for p in range(0,v,piv):
        pump.write(1)
        time.sleep(0.25)
        pump.write(0)
        time.sleep(0.25)
    pumpOnBtn.config(state=NORMAL)

def pumponvol():
    psGoBtn.config(state=DISABLED)
    cleanBtn.config(state=DISABLED)
    washBtn.config(state=DISABLED)
    resetBtn.config(state=DISABLED)
    vol = int(pumpvol.get())
    pumpon(vol)
    psGoBtn.config(state=NORMAL)
    cleanBtn.config(state=NORMAL)
    washBtn.config(state=NORMAL)
    resetBtn.config(state=NORMAL)

def go():
    pos = int(pspos.get())
    position = 'GO%d\r' % (pos)
    ps.write(position.encode())

def clean():
    pos1 = int(ps1pos.get())
    pos2 = int(ps2pos.get())
    psGoBtn.config(state=DISABLED)
    cleanBtn.config(state=DISABLED)
    washBtn.config(state=DISABLED)
    resetBtn.config(state=DISABLED)
    for n in range(pos1, pos2+1):
        cleaninfo = 'Cleaning port ' + str(n)
        cleanstatus.set(cleaninfo)
        position = 'GO%d\r' % (n)
        ps.write(position.encode())
        primeOn()
        pumpon(500+len1+len2)
        primeOff()
    ps.write('HM\r'.encode())
    cleanstatus.set('Cleaning completed')
    psGoBtn.config(state=NORMAL)
    cleanBtn.config(state=NORMAL)
    washBtn.config(state=NORMAL)
    resetBtn.config(state=NORMAL)

def wash():
    x = int(washTimes.get())
    psGoBtn.config(state=DISABLED)
    cleanBtn.config(state=DISABLED)
    washBtn.config(state=DISABLED)
    resetBtn.config(state=DISABLED)
    for n in range(1, x+1):
        washinfo = 'Resin washing ' + str(n)
        washstatus.set(washinfo)
        reagentOn()
        ps.write('GO2\r'.encode())
        pumpon(2000)
        ps.write('HM\r'.encode())
        time.sleep(1)
        reagentOff()
        n2On()
        wasteOn()
        ventOn()
        time.sleep(60)
        ventOff()
        wasteOff()
        n2Off()
    washstatus.set('Resin washing completed')
    psGoBtn.config(state=NORMAL)
    cleanBtn.config(state=NORMAL)
    washBtn.config(state=NORMAL)
    resetBtn.config(state=NORMAL)

def reset():
    volvar.set('0')
    psvar.set('1')
    ps1var.set('1')
    ps2var.set('1')
    washvar.set('0')
    cleanstatus.set('Ready')
    washstatus.set('Ready')
    n2Off()
    ventOff()
    reagentOff()
    wasteOff()
    primeOff()
    go()
    psGoBtn.config(state=NORMAL)
    cleanBtn.config(state=NORMAL)
    washBtn.config(state=NORMAL)
    resetBtn.config(state=NORMAL)
        
# Main
config = configparser.ConfigParser()
config.readfp(open('config.txt'))
pscom = config.get('Parameters', 'pscom')
arduinocom = config.get('Parameters', 'arduinocom')
tubevol = float(config.get('Parameters', 'tubevol'))
length1 = float(config.get('Parameters', 'length1'))
length2 = float(config.get('Parameters', 'length2'))
piv = int(config.get('Parameters', 'piv'))

ps = serial.Serial(port=pscom, baudrate=9600, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, bytesize=serial.EIGHTBITS) # VICI port selector
board = Arduino(arduinocom) # Arduino Uno
n2 = board.get_pin('d:2:o') # Solenoid valve normally closed, connected to digital pin 2
vent = board.get_pin('d:3:o') # Solenoid valve normally open, connected to digital pin 3
reagent = board.get_pin('d:4:o') # Solenoid valve normally closed, connected to digital pin 4
waste = board.get_pin('d:5:o') # Solenoid valve normally closed, connected to digital pin 5
prime = board.get_pin('d:6:o') # Solenoid valve normally closed, connected to digital pin 6
pump = board.get_pin('d:7:o') # Solenoid micro pump with an internal volume of 20 microliter and maximum pumping speed 2.4 ml/min or 40 microliter/sec, connected to digital pin 7 

font16 = ('Helvetica', 16, 'bold')
font12 = ('Helvetica', 12, 'bold')

root = Tk()
root.wm_title('PepSy')
frame = Frame(root)
frame.grid()

len1 = int(tubevol*length1) # tubing volume aa to ps
len2 = int(tubevol*length1) # tubing volume ps to pump

volvar = StringVar()
psvar = StringVar()
ps1var = StringVar()
ps2var = StringVar()
washvar = StringVar()
cleanstatus = StringVar()
washstatus = StringVar()

line1Lbl = Label(frame, text = '-----------------------------------------------------------------')
line1Lbl.grid(row = 1, column = 1, columnspan = 3)

title1Lbl = Label(frame, text = 'PepSy', font = font16)
title1Lbl.grid(row = 2, column = 1, columnspan = 3) 

line2Lbl = Label(frame, text = '-----------------------------------------------------------------')
line2Lbl.grid(row = 3, column = 1, columnspan = 3)

title2Lbl = Label(frame, text = 'Control Panel', font = font12)
title2Lbl.grid(row = 4, column = 1, columnspan = 3)

n2Lbl = Label(frame, text = 'Nitrogen Gas')
n2Lbl.grid(row = 5, column = 1, pady = 5, sticky = E)
n2OnBtn = Button(frame, text = ' ON ', command = lambda:n2On())
n2OnBtn.grid(row = 5, column = 2, pady = 5, sticky = E)
n2OffBtn = Button(frame, text = ' OFF ', command = lambda:n2Off())
n2OffBtn.grid(row = 5, column = 3, pady = 5, sticky = W)

ventLbl = Label(frame, text = 'Vent')
ventLbl.grid(row = 6, column = 1, pady = 5, sticky = E)
ventOnBtn = Button(frame, text = ' ON ', command = lambda:ventOn())
ventOnBtn.grid(row = 6, column = 2, pady = 5, sticky = E)
ventOffBtn = Button(frame, text = ' OFF ', command = lambda:ventOff())
ventOffBtn.grid(row = 6, column = 3, pady = 5, sticky = W)

reagentLbl = Label(frame, text = 'Reagent')
reagentLbl.grid(row = 7, column = 1, pady = 5, sticky = E)
reagentOnBtn = Button(frame, text = ' ON ', command = lambda:reagentOn())
reagentOnBtn.grid(row = 7, column = 2, pady = 5, sticky = E)
reagentOffBtn = Button(frame, text = ' OFF ',command = lambda:reagentOff())
reagentOffBtn.grid(row = 7, column = 3, pady = 5, sticky = W)

wasteLbl = Label(frame, text = 'Waste')
wasteLbl.grid(row = 8, column = 1, pady = 5, sticky = E)
wasteOnBtn = Button(frame, text = ' ON ', command = lambda:wasteOn())
wasteOnBtn.grid(row = 8, column = 2, pady = 5, sticky = E)
wasteOffBtn = Button(frame, text = ' OFF ', command = lambda:wasteOff())
wasteOffBtn.grid(row = 8, column = 3, pady = 5,sticky = W)

primeLbl = Label(frame, text = 'Prime')
primeLbl.grid(row = 9, column = 1, pady = 5, sticky = E)
primeOnBtn = Button(frame, text = ' ON ', command = lambda:primeOn())
primeOnBtn.grid(row = 9, column = 2, pady = 5, sticky = E)
primeOffBtn = Button(frame, text = ' OFF ',command = lambda:primeOff())
primeOffBtn.grid(row = 9, column = 3, pady = 5, sticky = W)

pumpLbl = Label(frame, text = 'Pump Volume')
pumpLbl.grid(row = 10, column = 1, pady = 5, sticky = E)
pumpvol = Entry(frame, textvariable = volvar, width = 5)
pumpvol.grid(row = 10, column = 2, pady = 5, sticky = E)
pumpOnBtn = Button(frame, text = ' ON ', command = lambda:_thread.start_new(pumponvol, ()))
pumpOnBtn.grid(row = 10, column = 3, pady = 5, sticky = W)

psLbl = Label(frame, text = 'Port')
psLbl.grid(row = 11, column = 1, pady = 5, sticky = E)
pspos = Entry(frame, textvariable = psvar, width = 5)
pspos.grid(row = 11, column = 2, pady = 5, sticky = E)
psGoBtn = Button(frame, text = ' GO ', command = lambda:go())
psGoBtn.grid(row = 11, column = 3, pady = 5, sticky = W)

line3Lbl = Label(frame, text = '-----------------------------------------------------------------')
line3Lbl.grid(row = 12, column = 1, columnspan = 3)

title4Lbl = Label(frame, text = 'Amino Acid/Reagent Line Cleaning', font = font12)
title4Lbl.grid(row = 13, column = 1, columnspan = 3)

title5Lbl = Label(frame, text = ' ')
title5Lbl.grid(row = 14, column = 1, columnspan = 3)

ps1Lbl = Label(frame, text = 'Starting Port')
ps1Lbl.grid(row = 15, column = 1, sticky = E)
ps1pos = Entry(frame, textvariable = ps1var, width = 5)
ps1pos.grid(row = 15, column = 2, sticky = E)

ps2Lbl = Label(frame, text = 'Ending Port')
ps2Lbl.grid(row = 16, column = 1, sticky = E)
ps2pos = Entry(frame, textvariable = ps2var, width = 5)
ps2pos.grid(row = 16, column = 2, sticky = E)

cleanBtn = Button(frame, text = '          CLEAN          ', command = lambda:_thread.start_new(clean, ()))
cleanBtn.grid(row = 17, column = 1, pady = 5, columnspan = 3)

cleanstatusLbl = Label(frame, textvariable = cleanstatus)
cleanstatusLbl.grid(row = 18, column = 1, pady = 5, columnspan = 3)

line4Lbl = Label(frame, text = '-----------------------------------------------------------------')
line4Lbl.grid(row = 19, column = 1, columnspan = 3)

title6Lbl = Label(frame, text = 'Resin Washing', font = font12)
title6Lbl.grid(row = 20, column = 1, columnspan = 3)

title7Lbl = Label(frame, text = ' ')
title7Lbl.grid(row = 21, column = 1, columnspan = 3)

washLbl = Label(frame, text = 'Times')
washLbl.grid(row = 22, column = 1, sticky = E)
washTimes = Entry(frame, textvariable = washvar, width = 5)
washTimes.grid(row = 22, column = 2, sticky = E)

washBtn = Button(frame, text = '          WASH          ', command = lambda:_thread.start_new(wash, ()))
washBtn.grid(row = 23, column = 1, pady = 5, columnspan = 3)

washstatusLbl = Label(frame, textvariable = washstatus)
washstatusLbl.grid(row = 24, column = 1, pady = 5, columnspan = 3)

line5Lbl = Label(frame, text = '-----------------------------------------------------------------')
line5Lbl.grid(row = 25, column = 1, columnspan = 3)

resetBtn = Button(frame, text = '       RESET          ', command = lambda:reset())
resetBtn.grid(row = 26, column = 1, pady = 5, columnspan = 3)

reset()

root.mainloop()
