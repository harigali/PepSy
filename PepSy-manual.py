# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PepSy - An open-source peptide synthesizer
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Developed by Dr. Hariprasad Gali, Ph.D., Research Assistant Professor, Department of Pharmaceutical Sciences, College of Pharmacy, The University of Oklahoma Health Sciences Center, Oklahoma City, OK 73117.
# Email address to report bugs: hgali@ouhsc.edu.

# Tested only with Python 3.5.0
# First full version - December 16, 2015
# Last update - May 12, 2016

# This script is written for the manual control of PepSy and cleaning amino acid/reagent line.
# Save device configuration file (config.txt) in the same folder where PepSy.py and PepSy-manual.py scripts are saved.
# Only Arduino UNO digital pins are used.
# COM port numbers of VICI stream selector valve (ps) and Arduino UNO (board) needs to be updated in the device configuration file according to their current assignment on the PC.

# VICI CHEMINERT low pressure 24 position stream selector valve
# Position 1 - Air
# Position 2 - DMF
# Position 3 - DCM
# Position 4 - Piperidine solution
# Position 5 - DIPEA solution
# Position 6 - HOBT solution
# Position 7 - HBTU solution
# Position 8 to 24 - Amino acid solutions or other reagent solutions

import time
import serial
import configparser
from tkinter import *
from pyfirmata import Arduino, util

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

root = Tk()
root.wm_title("PepSy")
frame = Frame(root)
frame.grid()

len1 = int(tubevol*length1) # tubing volume aa to ps
len2 = int(tubevol*length1) # tubing volume ps to pump

# Functions
def pumpon(v): # v is volume (integer) to be pumped in microliters, pumping rate is 1.2 ml/min, this function is required as a solenoid valve based micro pump is being used
    for p in range(0,v,piv):
        pump.write(1)
        time.sleep(0.25)
        pump.write(0)
        time.sleep(0.25)

def on():
    vol = int(pin7vol.get())
    pumpon(vol)

def go():
    pos = int(pspos.get())
    position = 'GO%d\r' % (pos)
    ps.write(position.encode())

def clean():
    pos1 = int(ps1pos.get())
    pos2 = int(ps2pos.get())
    for n in range(pos1, pos2+1):
        position = 'GO%d\r' % (n)
        ps.write(position.encode())
        prime.write(1)
        pumpon(500+len1+len2)
        prime.write(0)
    ps.write("HM\r".encode())
    status.set("Cleaning completed")

def washing():
    x = int(washTimes.get())
    for n in range(1, x+1):
        reagent.write(1)
        ps.write("GO2\r".encode())
        pumpon(2000)
        ps.write("HM\r".encode())
        time.sleep(1)
        reagent.write(0)
        n2.write(1)
        waste.write(1)
        vent.write(1)
        time.sleep(60)
        vent.write(0)
        waste.write(0)
        n2.write(0)
    washstatus.set("Resin washing completed")

# Main
line1Lbl = Label(frame, text = "-----------------------------------------------------------------")
line1Lbl.grid(row = 1, column = 1, columnspan = 3)

title1Lbl = Label(frame, text = "PepSy")
title1Lbl.grid(row = 2, column = 1, columnspan = 3)

line2Lbl = Label(frame, text = "-----------------------------------------------------------------")
line2Lbl.grid(row = 3, column = 1, columnspan = 3)

title2Lbl = Label(frame, text = "Control Panel")
title2Lbl.grid(row = 4, column = 1, columnspan = 3)

pin2Lbl = Label(frame, text = "Nitrogen Gas")
pin2Lbl.grid(row = 5, column = 1, pady = 5, sticky = E)
pin2onbutton = Button(frame, text = " ON ", command = lambda:n2.write(1))
pin2onbutton.grid(row = 5, column = 2, pady = 5, sticky = E)
pin2offbutton = Button(frame, text = " OFF ", command = lambda:n2.write(0))
pin2offbutton.grid(row = 5, column = 3, pady = 5, sticky = W)

pin3Lbl = Label(frame, text = "Vent")
pin3Lbl.grid(row = 6, column = 1, pady = 5, sticky = E)
pin3onbutton = Button(frame, text = " ON ", command = lambda:vent.write(1))
pin3onbutton.grid(row = 6, column = 2, pady = 5, sticky = E)
pin3offbutton = Button(frame, text = " OFF ", command = lambda:vent.write(0))
pin3offbutton.grid(row = 6, column = 3, pady = 5, sticky = W)

pin4Lbl = Label(frame, text = "Reagent")
pin4Lbl.grid(row = 7, column = 1, pady = 5, sticky = E)
pin4onbutton = Button(frame, text = " ON ", command = lambda:reagent.write(1))
pin4onbutton.grid(row = 7, column = 2, pady = 5, sticky = E)
pin4offbutton = Button(frame, text = " OFF ",command = lambda:reagent.write(0))
pin4offbutton.grid(row = 7, column = 3, pady = 5, sticky = W)

pin5Lbl = Label(frame, text = "Waste")
pin5Lbl.grid(row = 8, column = 1, pady = 5, sticky = E)
pin5onbutton = Button(frame, text = " ON ", command = lambda:waste.write(1))
pin5onbutton.grid(row = 8, column = 2, pady = 5, sticky = E)
pin5offbutton = Button(frame, text = " OFF ", command = lambda:waste.write(0))
pin5offbutton.grid(row = 8, column = 3, pady = 5,sticky = W)

pin6Lbl = Label(frame, text = "Prime")
pin6Lbl.grid(row = 9, column = 1, pady = 5, sticky = E)
pin6onbutton = Button(frame, text = " ON ", command = lambda:prime.write(1))
pin6onbutton.grid(row = 9, column = 2, pady = 5, sticky = E)
pin6offbutton = Button(frame, text = " OFF ",command = lambda:prime.write(0))
pin6offbutton.grid(row = 9, column = 3, pady = 5, sticky = W)

pin7Lbl = Label(frame, text = "Pump Volume")
pin7Lbl.grid(row = 10, column = 1, pady = 5, sticky = E)
pin7vol = Entry(frame, width = 5)
pin7vol.grid(row = 10, column = 2, pady = 5, sticky = E)
pin7onbutton = Button(frame, text = " ON ", command = lambda:on())
pin7onbutton.grid(row = 10, column = 3, pady = 5, sticky = W)

psLbl = Label(frame, text = "Port")
psLbl.grid(row = 11, column = 1, pady = 5, sticky = E)
pspos = Entry(frame, width = 5)
pspos.grid(row = 11, column = 2, pady = 5, sticky = E)
psgobutton = Button(frame, text = " GO ", command = lambda:go())
psgobutton.grid(row = 11, column = 3, pady = 5, sticky = W)

line3Lbl = Label(frame, text = "-----------------------------------------------------------------")
line3Lbl.grid(row = 12, column = 1, columnspan = 3)

title4Lbl = Label(frame, text = "Amino acid/reagent Line Cleaning")
title4Lbl.grid(row = 13, column = 1, columnspan = 3)

title5Lbl = Label(frame, text = " ")
title5Lbl.grid(row = 14, column = 1, columnspan = 3)

ps1Lbl = Label(frame, text = "Starting Port")
ps1Lbl.grid(row = 15, column = 1, sticky = E)
ps1pos = Entry(frame, width = 5)
ps1pos.grid(row = 15, column = 2, sticky = E)

ps2Lbl = Label(frame, text = "Ending Port")
ps2Lbl.grid(row = 16, column = 1, sticky = E)
ps2pos = Entry(frame, width = 5)
ps2pos.grid(row = 16, column = 2, sticky = E)

ps1gobutton = Button(frame, text = "          CLEAN          ", command = lambda:clean())
ps1gobutton.grid(row = 17, column = 1, pady = 5, columnspan = 3)

status = StringVar()
statusLbl = Label(frame, textvariable = status)
statusLbl.grid(row = 18, column = 1, pady = 5, columnspan = 3)

line4Lbl = Label(frame, text = "-----------------------------------------------------------------")
line4Lbl.grid(row = 19, column = 1, columnspan = 3)

title6Lbl = Label(frame, text = "Resin Washing")
title6Lbl.grid(row = 20, column = 1, columnspan = 3)

title7Lbl = Label(frame, text = " ")
title7Lbl.grid(row = 21, column = 1, columnspan = 3)

washLbl = Label(frame, text = "Times")
washLbl.grid(row = 22, column = 1, sticky = E)
washTimes = Entry(frame, width = 5)
washTimes.grid(row = 22, column = 2, sticky = E)

ps1gobutton = Button(frame, text = "          Wash          ", command = lambda:washing())
ps1gobutton.grid(row = 23, column = 1, pady = 5, columnspan = 3)

washstatus = StringVar()
washstatusLbl = Label(frame, textvariable = washstatus)
washstatusLbl.grid(row = 24, column = 1, pady = 5, columnspan = 3)

line5Lbl = Label(frame, text = "-----------------------------------------------------------------")
line5Lbl.grid(row = 25, column = 1, columnspan = 3)

root.mainloop()
