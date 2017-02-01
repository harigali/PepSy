# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PepSy - An open-source peptide synthesizer
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Developed by Dr. Hariprasad Gali, Ph.D., Research Assistant Professor, Department of Pharmaceutical Sciences, College of Pharmacy, The University of Oklahoma Health Sciences Center, Oklahoma City, OK 73117.
# Email address to report bugs: hgali@ouhsc.edu.
# Tested only with Python 3.5.0
# Last update - June 8, 2016

# Create folders named "sequence" and "output" within the same folder where PepSy.py and PepSy-manual.py scripts are saved.
# Save device configuration file (config.txt) in the same folder where PepSy.py and PepSy-manual.py scripts are saved.
# Create a sequence configuration file (see example templete.txt) for each run and save it in the "sequence" folder.
# An output file is generated for each run and stored in the "output" folder.

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

# Amino acid priming requires 0.30 ml for single coupling and 0.43 ml for double coupling.
# Piperidine, DIPEA, HOBT, and HBTU priming requires 0.5 ml initially and 0.3 ml for each amino acid and other reagents.
# Amino acid priming volumes are optimized for tube lengths (~11.6 microliter per inch for 0.03" ID tubing) - ~15 inch from aa to ps, ~11 inch from ps to pump, and ~20 inch pump to resin.
# Tubing volume variables are len1 (aa to ps),  len2 (ps to pump), and  len3 (pump to resin) depds on tubing id and lengths, which needs to be updated in the device configuration file when replacing the tubing with different length or id.
# Solution concentrations to be used: amino acids - 0.33M, HBTU - 0.33M, HOBT - 0.66M, piperidine - 20%, hydrazine - 2%, acetic anhydride/pyridine - 2.5M/2.5M, and Tl(CF3COO)2 - 0.05M.
# If the solenoid micro pump is replaced with a different pump internal volume (piv), it needs to be updated in the device configuration file.

# Uppercase alphabets are used for both L and D amino acids.
# Lowercase alphabets are used for N-methyl amino acids.
# "X", "1", "2", "3", "4", or "5" are used for a linker, an unusual amino acid or any molecule that requires both coupling and fmoc deprotection and place the solution in the position assigned to "X", "1", "2", "3", "4", or "5" respectively.
# "Z", "6", "7", "8", or "9" are used for a chelator, an unusual amino acid, or any molecule that requires only coupling and place the solution in the position assigned to "Z", "6", "7", "8", or "9" respectively.
# "*" is used for pausing the synthesis.
# "!" is used for ivDde deprotection and place the hydrazine solution in the position assigned to "!".
# "@" is used for onresin oxidation and place the thallium solution in the position assigned to "@".
# "$" is used for endcapping and place the acetic anhydride solution in the position assigned to "$".
# -------------------------------------------------------------------------------------------------------------------------------------------

# Imports
import os
import time
import datetime
import serial
import configparser
from pyfirmata import Arduino, util
from collections import Counter
# -------------------------------------------------------------------------------------------------------------------------------------------

# Functions
def positions(p):
    mwdict = {"A":329.36, "a":329.36, "C":585.72, "c":585.72, "D":411.45, "d":411.45, "E":425.48, "e":425.48, "F":387.44, "f":387.44, "G":297.31, "g":297.31, "H":619.72, "h":619.72, "I":353.42, "i":353.42, "K":468.2, "k":468.2,
				"L":353.42, "l":353.42, "M":371.45, "m":371.45, "N":596.68, "n":596.68, "P":337.38, "p":337.38, "Q":610.71, "q":610.71, "R":648.78, "r":648.78, "S":383.44, "s":383.44, "T":379.48, "t":379.48, "V":339.39, "v":339.39,
				"W":526.59, "w":526.59, "Y":459.54, "y":459.54} # molecular weight of standard fmoc-protected amino acids
    paap = [] # positions for different amino acids and reagents
    pseq = Counter(x for x in p if x not in ignore) # amino acids and reagents sorting
    paan1 = len(pseq) # number of different amino acids and reagents
    paak = list(pseq.keys())
    paav = list(pseq.values())
    file = open(filename, 'a')
    if pa == "y" or pa == "Y":
        print("Place amino acid/reagent solutions with required volumes in the positions shown below")
        print(" ")
        print("---------------------------------------------------------------------------------------------------")
        file.write("---------------------------------------------------------------------------------------------------" + '\n')
        print("S. No.", '\t', "Amino acid", '\t', "Position", '\t', "Solution volume", '\t', "Amino acid weight", '\t', "DMF volume")
        print("---------------------------------------------------------------------------------------------------")
        file.write("S. No." + '\t' + "Amino acid" + '\t' + "Position" + '\t' + "Solution volume" + '\t\t' + "Amino acid weight" + '\t' + "DMF volume" + '\n')
        file.write("---------------------------------------------------------------------------------------------------" + '\n')
        for n in range (1, paan1+1):
            at = n % (ports - 7)
            if at == 0:
                pos = ports
            else:
                pos = 7 + at
            paap.append(pos)
            if paak[n-1] == "P" or paak[n-1].islower():
                vol = ((len1+len2*2)/1000+ss*1)*paav[n-1]
            elif paak[n-1] == "!":
                vol = (len1+len2)/1000+ss*2
            elif paak[n-1] == "@":
                vol = (len1+len2)/1000+ss*4
            elif paak[n-1] == "$":
                vol = (len1+len2)/1000+ss*1
            else:
                vol = ((len1+len2)/1000+ss*0.5)*paav[n-1]  
            try:
                mw = mwdict[paak[n-1]]
            except:
                mw = 0
            wt = vol*mw*0.33
            dmf = vol*1000-wt
            vol = str(round(vol, 2))
            wt = format(int(wt), '03d')
            dmf = str(int(dmf))
            print(n, '\t', paak[n-1], "(", paav[n-1], ")", '\t', paap[n-1], '\t\t', vol, "ml", '\t\t', wt, "mg", '\t\t', dmf, "ul") # displays wt for standard amino acids
            file.write(str(n) + '\t' + paak[n-1] + "(" + str(paav[n-1]) + ")" + '\t\t' + str(paap[n-1]) + '\t\t' + vol + " ml" + '\t\t\t' + wt + " mg" + '\t\t\t' + dmf + " ul" + '\n')
        print("---------------------------------------------------------------------------------------------------")
        file.write("---------------------------------------------------------------------------------------------------" + '\n')
        file.close()
        print(" ")
        print("Check the nitrogen gas pressure, if it is not ~2 psi then adjust the pressure.")
        print(" ")
        print("Check the levels of all the reagents, if any of them is not enough then add.")
        print(" ")
        input("If you are ready, press ENTER to continue")
    else:
        for n in range (1, paan1+1):
            pos = synconfig.getint('Positions', paak[n-1])
            paap.append(pos)
    print(" ")
    print("----------------------------------------------------------------------")
    file = open(filename, 'a')
    file.write("----------------------------------------------------------------------" + '\n')
    print("S. No.", '\t', "Amino acid", '\t', "Position", '\t', "Coupling", '\t', "Deprotection")
    print("----------------------------------------------------------------------")
    file.write("S. No." + '\t' + "Amino acid" + '\t' + "Position" + '\t' + "Coupling" + '\t' + "Deprotection" + '\n')
    file.write("----------------------------------------------------------------------" + '\n')
    for n in range (1, paan+1):
        if aa[n-1] == "*":
            a.append(1)
        for m in range (1, paan1+1):
            if aa[n-1] == paak[m-1]:
                a.append(paap[m-1])
        if aa[n-1] == "P" or aa[n-1].islower():
            c.append("double") # double coupling if aa is P or any aa represented by a lowercase letter
        elif aa[n-1] == "!":
            c.append("ivdde")
        elif aa[n-1] == "@":
            c.append("oxidation")
        elif aa[n-1] == "$":
            c.append("encapping")
        elif aa[n-1] == "*":
            c.append("pause")
        else:
            c.append("single") # default coupling
        if aa[n-1] == "*" or aa[n-1] == "!" or aa[n-1] == "@" or aa[n-1] == "$" or aa[n-1] == "Z" or aa[n-1] == "6" or aa[n-1] == "7" or aa[n-1] == "8" or aa[n-1] == "9":
            d.append("none")
        else:
            d.append("fmoc") # default deprotection
        if saa > 1:
            n1 = n + saa - 1
        else:
            n1 = n
        print(n1, '\t', aa[n-1], '\t\t', a[n-1], '\t\t', c[n-1], '\t', d[n-1])
        file.write(str(n1) + '\t' + aa[n-1] + '\t\t' + str(a[n-1]) + '\t\t' + c[n-1] + '\t\t' + d[n-1] + '\n')
    print("----------------------------------------------------------------------")
    print(" ")
    file.write("----------------------------------------------------------------------" + '\n')
    file.close()
    input("If the positions, couplings, and deprotections are correct, press ENTER to continue")
    print(" ")
    
def pumpon(v): # v is volume (integer) to be pumped in microliters
    for p in range(0, v, piv): 
        pump.write(1)
        time.sleep(0.25)
        pump.write(0)
        time.sleep(0.25)

def presyn():
    initialization()
    if pr == "y" or pr == "Y":
        priming()
    elif pr == "n" or pr == "N":
        print("Priming skipped")
        print(" ")
    else:
        pr1 = input("Input error in the sequence file. Do you want to perform priming (y or n)? ")
        print(" ")
        if pr1 == "y" or pr1 == "Y":
            priming()
        else:
            print("Priming skipped")
            print(" ")
    if sw == "y" or sw == "Y":
        swelling()
    elif sw == "n" or sw == "N":
        print("Swelling skipped")
        print(" ")
    else:
        sw1 = input("Input error in the sequence file. Do you want to perform swelling (y or n)? ")
        print(" ")
        if sw1 == "y" or sw1 == "Y":
            swelling()
        else:
            print("Swelling skipped")
            print(" ")
    print("--------------------------------------------------------")
    print("Synthesis")
    print("--------------------------------------------------------")
    file = open(filename, 'a')
    file.write("Peptide synthesis started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    if dp == "y" or dp == "Y":
        fmocdeprotection()
    elif dp == "n" or dp == "N":
        print("Initial fmoc deprotection skipped")
        print(" ")
    else:
        dp1 = input("Input error in the sequence file. Do you want to perform initial fmoc deprotection (y or n)? ")
        print(" ")
        if dp1 == "y" or dp1 == "Y":
            fmocdeprotection()
        else:
            print("Initial fmoc deprotection skipped")
            print(" ")
        
def syn():
    for n in range(1, paan+1):
        if saa > 1:
            n1 = n + saa - 1
        else:
            n1 = n
        if aa[n-1] == "*":
            print("Pause")
            print("------------------")
            file = open(filename, 'a')
            file.write("Pause" + '\n')
            file.close()
        elif aa[n-1] == "@":
            print("On-resin oxidation")
            print("------------------")
            file = open(filename, 'a')
            file.write("On-resin oxidation" + '\n')
            file.close()
        elif aa[n-1] == "!":
            print("ivDde deprotection")
            print("------------------")
            file = open(filename, 'a')
            file.write("ivDde deprotection" + '\n')
            file.close()
        elif aa[n-1] == "$":
            print("Endcapping")
            print("------------------")
            file = open(filename, 'a')
            file.write("Endcapping" + '\n')
            file.close()
        else:
            print("Amino acid: " + str(n1) + " (" + aa[n-1] + ")")
            print("------------------")
            file = open(filename, 'a')
            file.write("Amino acid: " + str(n1) + " (" + aa[n-1] + ")" + '\n')
            file.close()
        if c[n-1] == "single": 
            coupling(n-1)
        elif c[n-1] == "double":
            doublecoupling(n-1)
        elif c[n-1] == "oxidation":
            onresinoxidation(n-1)
        elif c[n-1] == "endcapping":
            endcapping(n-1) 
        elif c[n-1] == "ivdde":
            ivddedeprotection(n-1) 
        elif c[n-1] == "pause":
            pause()
        else:
            print("No coupling")
            print(" ")
        if d[n-1] == "fmoc":
            fmocdeprotection()
    
def initialization():
    print("--------------------------------------------------------")
    print("Initialization")
    print("--------------------------------------------------------")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(0)
    vent.write(0)
    reagent.write(0)
    waste.write(0)
    prime.write(0)
    pump.write(0)
    file = open(filename, 'a')
    file.write("Initialization completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def priming():
    print("--------------------------------------------------------")
    print("Priming")
    print("--------------------------------------------------------")
    print("DMF, DCM, piperidine, DIPEA, HOBT, and HBTU lines")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    prime.write(1)
    for p in range(4,8):
        position = 'GO%d\r' % (p)
        ps.open()
        ps.write(position.encode())
        ps.close()
        pumpon(len1+len2)
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        time.sleep(1)
    for p in range(2,4):
        position = 'GO%d\r' % (p)
        ps.open()
        ps.write(position.encode())
        ps.close()
        pumpon(1000)
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        time.sleep(1)
    prime.write(0)
    file = open(filename, 'a')
    file.write("Priming completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def swelling():
    print("--------------------------------------------------------")
    print("Swelling")
    print("--------------------------------------------------------")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print("Adding solvents")
    reagent.write(1)
    ps.open()
    ps.write("GO3\r".encode())
    ps.close()
    pumpon(1000-len3) # addition of 1 ml DCM, 768 from pumping, 232 from tubing
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    time.sleep(1)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(1000+len3) # addition of 1 ml DMF, 1000 DMF, 232 residual DCM
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    n2.write(1)
    print("15 min swelling")
    time.sleep(900) # 15 min swelling
    print("Draining solvents")
    waste.write(1)
    vent.write(1)
    time.sleep(30) # draining
    n2.write(0)
    waste.write(0)
    vent.write(0)
    print("Washing")
    washing() # washing
    file = open(filename, 'a')
    file.write("Swelling completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")
    
def coupling(n): # S. No. (integer) of the amino acid
    file = open(filename, 'a')
    file.write("Coupling (single) started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    print("Coupling (single)")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    aapos = a[n]
    position = 'GO%d\r' % (aapos)
    ps.open()
    ps.write(position.encode())
    ps.close()
    prime.write(1)
    file.write("Amino acid position on PS is " + str(aapos) + '\n')
    file.close()
    print("Amino acid position on PS is " + str(aapos))
    print("Priming amino acid " + aa[n])
    pumpon(len1+len2) # amino acid line priming - aa to ps to pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # amino acid line priming - pump to resin
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    print("Adding reagents")
    reagent.write(1)
    pumpon(ss*500-len3) # addition of 0.5 ml amino acid solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    time.sleep(1)
    ps.open()
    ps.write("GO5\r".encode())
    ps.close()
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(ss*260) # addition of 0.26 ml DIPEA solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    time.sleep(1)
    ps.open()
    ps.write("GO6\r".encode())
    ps.close()
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(ss*260) # addition of 0.26 ml HOBT solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    time.sleep(1)
    ps.open()
    ps.write("GO7\r".encode())
    ps.close()
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(ss*500) # addition of 0.5 ml HBTU solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing pump to resin
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    n2.write(1)
    print("60 min coupling")
    time.sleep(3600) # 60 min coupling
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    n2.write(0)
    waste.write(0)
    vent.write(0)
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        washing()
    file = open(filename, 'a')
    file.write("Coupling completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def doublecoupling(n): # S. No. (integer) of the amino acid
    file = open(filename, 'a')
    file.write("Coupling (double) started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Coupling (double)")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    for d in range(1,3):
        aapos = a[n]
        position = 'GO%d\r' % (aapos)
        ps.open()
        ps.write(position.encode())
        ps.close()
        prime.write(1)
        if d == 1:
            file = open(filename, 'a')
            file.write("Amino acid position on PS is " + str(aapos) + '\n')
            file.close()
            print("Amino acid position on PS is " + str(aapos))
            print("Priming amino acid " + aa[n])
            pumpon(len1+len2) # amino acid line priming - aa to ps to pump
        if d == 2:
            pumpon(len2) # amino acid line priming - ps to pump
        prime.write(0)
        reagent.write(1)
        pumpon(len3) # amino acid line priming - pump to resin
        reagent.write(0)
        waste.write(1)
        vent.write(1)
        n2.write(1)
        time.sleep(10)
        n2.write(0)
        vent.write(0)
        waste.write(0)
        print("Adding reagents")
        reagent.write(1)
        pumpon(ss*500-len3) # addition of 0.5 ml amino acid solution
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        time.sleep(1)
        ps.open()
        ps.write("GO5\r".encode())
        ps.close()
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # removing previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(ss*260) # addition of 0.26 ml DIPEA solution
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        time.sleep(1)
        ps.open()
        ps.write("GO6\r".encode())
        ps.close()
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # removing previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(ss*260) # addition of 0.26 ml HOBT solution
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        time.sleep(1)
        ps.open()
        ps.write("GO7\r".encode())
        ps.close()
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # removing previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(ss*500) # addition of 0.5 ml HBTU solution
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        reagent.write(0)
        ps.open()
        ps.write("GO2\r".encode())
        ps.close()
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(len3) # DMF to add previous reagent leftover in the tubing
        reagent.write(0)
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        n2.write(1)
        print("60 min coupling")
        time.sleep(3600) # 60 min coupling
        waste.write(1)
        vent.write(1)
        print("Draining reagents")
        time.sleep(30)# draining
        n2.write(0)
        waste.write(0)
        vent.write(0)
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        washing()
    file = open(filename, 'a')
    file.write("Double coupling completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def fmocdeprotection():
    file = open(filename, 'a')
    file.write("fmoc deprotection started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("fmoc deprotection")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print("Adding reagents")
    ps.open()
    ps.write("GO4\r".encode())
    ps.close()
    time.sleep(1)
    prime.write(1)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*1000) # addition of 1 ml piperidine solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("10 min first round deprotection")
    time.sleep(600) # 10 min first round deprotection
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    print("Adding reagents")
    reagent.write(1)
    ps.open()
    ps.write("GO4\r".encode())
    ps.close()
    pumpon(ss*1000) # addition of 1 ml piperidine solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("20 min second round deprotection")
    time.sleep(1200) # 20 min second round deprotection
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        washing()
    file = open(filename, 'a')
    file.write("Deprotection completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")
    
def washing():
    reagent.write(1)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    pumpon(2000) # addition of 2 ml DMF
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    time.sleep(1)
    reagent.write(0)
    n2.write(1)
    waste.write(1)
    vent.write(1)
    time.sleep(60) # draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    
def finalwashing():
    print("--------------------------------------------------------")
    print("Final washing")
    print("--------------------------------------------------------")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        reagent.write(1)
        ps.open()
        ps.write("GO3\r".encode())
        ps.close()
        pumpon(2000) # addition of 2 ml DCM
        ps.open()
        ps.write("HM\r".encode())
        ps.close()
        time.sleep(1)
        reagent.write(0)
        n2.write(1)
        waste.write(1)
        vent.write(1)
        time.sleep(60) # draining
        vent.write(0)
        waste.write(0)
        n2.write(0)
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")
    print("--------------------------------------------------------")
    print("Drying")
    print("--------------------------------------------------------")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    drying()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def pause():
    reagent.write(1)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    pumpon(1000) # addition of 1 ml DMF
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    input("Synthesis paused, Press ENTER to continue")
    print(" ")
    n2.write(1)
    waste.write(1)
    vent.write(1)
    time.sleep(15) # draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    
def aalinecleaning():
    m = int(input("Enter the starting position on ps "))
    print(" ")
    n = int(input("Enter the ending position on ps "))
    print(" ")
    input("Insert all amino acid/reagent lines in DMF and then press ENTER to continue")
    print(" ")
    for o in range(m,n+1):
        print("Cleaning line ", o)
        position = 'GO%d\r' % (o)
        ps.open()
        ps.write(position.encode())
        ps.close()
        prime.write(1)
        pumpon(500+len1+len2)
        prime.write(0)
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    print(" ")
    print("Remove amino acid/reagent lines from DMF and clean the exterior with acetone or isopropyl alcohol wipe")
    print(" ")
    print("Amino acid/reagent lines cleaning completed")
    print(" ")
    
def drying():
    n2.write(1)
    vent.write(1)
    waste.write(1)
    time.sleep(1800) # 30 min drying
    n2.write(0)
    vent.write(0)
    waste.write(0)

def ivddedeprotection(n): # S. No. (integer) of the reagent
    file = open(filename, 'a')
    file.write("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print("Adding reagents")
    aapos = a[n]
    position = 'GO%d\r' % (aapos)
    ps.open()
    ps.write(position.encode())
    ps.close()
    file.write("ivDde position on PS is " + str(aapos) + '\n')
    file.close()
    print("ivDde position on PS is " + str(aapos))
    time.sleep(1)
    prime.write(1)
    pumpon(len1+len2) # removing previous reagent from tubing between aa to ps to pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*1000) # addition of 1 ml hydrazine solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("10 min first round deprotection")
    time.sleep(600) # 10 min first round deprotection
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    print("Adding reagents")
    reagent.write(1)
    ps.open()
    ps.write(position.encode())
    ps.close()
    pumpon(ss*1000) # addition of 1 ml hydrazine solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("20 min second round deprotection")
    time.sleep(1200) # 20 min second round deprotection
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        washing()
    file = open(filename, 'a')
    file.write("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def onresinoxidation(n): # S. No. (integer) of the reagent
    file = open(filename, 'a')
    file.write("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print("Adding reagents")
    aapos = a[n]
    position = 'GO%d\r' % (aapos)
    ps.open()
    ps.write(position.encode())
    ps.close()
    file.write("Tl(CF3COO)3 position on PS is " + str(aapos) + '\n')
    file.close()
    print("Tl(CF3COO)3 position on PS is " + str(aapos))
    time.sleep(1)
    prime.write(1)
    pumpon(len1+len2) # removing previous reagent from tubing between aa to ps to pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*2000) # addition of 2 ml Tl(CF3COO)2 solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("60 min first round oxidation")
    time.sleep(3600) # 60 min first round oxidation
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    print("Adding reagents")
    reagent.write(1)
    ps.open()
    ps.write(position.encode())
    ps.close()
    pumpon(ss*2000) # addition of 2 ml Tl(CF3COO)2 solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("60 min second round oxidation")
    time.sleep(3600) # 60 min second round oxidation
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len2) # DMF to add previous reagent leftover in the tubing
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        washing()
    file = open(filename, 'a')
    file.write("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")

def endcapping(n): # S. No. (integer) of the reagent
    file = open(filename, 'a')
    file.write("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print("Adding reagents")
    aapos = a[n]
    position = 'GO%d\r' % (aapos)
    ps.open()
    ps.write(position.encode())
    ps.close()
    file.write("Acetic anhydride position on PS is " + str(aapos) + '\n')
    file.close()
    print("Acetic anhydride position on PS is " + str(aapos))
    time.sleep(1)
    prime.write(1)
    pumpon(len1+len2) # removing previous reagent from tubing between aa to ps to pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*1000) # addition of 1 ml acetic anhydride solution
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    n2.write(1)
    reagent.write(0)
    print("30 min end capping")
    time.sleep(1800) # 30 min end capping
    waste.write(1)
    vent.write(1)
    print("Draining reagents")
    time.sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    ps.open()
    ps.write("GO2\r".encode())
    ps.close()
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing
    ps.open()
    ps.write("HM\r".encode())
    ps.close()
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    time.sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print("Washing " + str(w))
        washing()
    file = open(filename, 'a')
    file.write("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p') + '\n')
    file.close()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    print(" ")
# -------------------------------------------------------------------------------------------------------------------------------------------

# Main
devconfig = configparser.ConfigParser()
devconfig.readfp(open('config.txt'))
pscom = devconfig.get('Parameters', 'pscom')
arduinocom = devconfig.get('Parameters', 'arduinocom')
ports = devconfig.getint('Parameters', 'ports')
tubevol = devconfig.getfloat('Parameters', 'tubevol')
length1 = devconfig.getfloat('Parameters', 'length1')
length2 = devconfig.getfloat('Parameters', 'length2')
length3 = devconfig.getfloat('Parameters', 'length3')
piv = devconfig.getint('Parameters', 'piv')

len1 = int(tubevol*length1) 
len2 = int(tubevol*length2) 
len3 = int(tubevol*length3)

ps = serial.Serial(port=pscom, baudrate=9600, parity=serial.PARITY_NONE, stopbits=serial.STOPBITS_ONE, bytesize=serial.EIGHTBITS) # VICI port selector
board = Arduino(arduinocom) # Arduino Uno
n2 = board.get_pin('d:2:o') # Solenoid valve normally closed, connected to digital pin 2
vent = board.get_pin('d:3:o') # Solenoid valve normally open, connected to digital pin 3
reagent = board.get_pin('d:4:o') # Solenoid valve normally closed, connected to digital pin 4
waste = board.get_pin('d:5:o') # Solenoid valve normally closed, connected to digital pin 5
prime = board.get_pin('d:6:o') # Solenoid valve normally closed, connected to digital pin 6
pump = board.get_pin('d:7:o') # Solenoid micro pump with an internal volume of 20 microliter and rated for a maximum pumping rate of 2.4 ml/min or 40 microliter/sec, connected to digital pin 7

print(" ")
print("--------------------------------------------------------------------------------------------------")
print("                                           PepSy                                                  ")
print("--------------------------------------------------------------------------------------------------")
print(" ")
print(datetime.datetime.now().strftime('%m-%d-%Y %I:%M:%S %p'))
print(" ")

seqfile = input("Enter the sequence configuration file name ")
print(" ")
if seqfile == "":
    q1 = input("Do you want to terminate the run (y or n)? ")
    print(" ")
    if q1 == "y" or q1 == "Y":
        ps.close()
        exit()
    else:
        seqfile = input("Enter the sequence configuration file name ")
        print(" ")
filename = seqfile + datetime.datetime.now().strftime('-%Y-%m-%d-') +  "out.txt"
seqfile = "sequence/" + seqfile + ".txt"

synconfig = configparser.ConfigParser()
synconfig.readfp(open(seqfile))
ss = synconfig.getint('Parameters', 'ss')
seq = synconfig.get('Parameters', 'seq')
pa = synconfig.get('Parameters', 'pa')
saa = synconfig.getint('Parameters', 'saa')
pr = synconfig.get('Parameters', 'pr')
sw = synconfig.get('Parameters', 'sw')
dp = synconfig.get('Parameters', 'dp')
fw = synconfig.get('Parameters', 'fw')
if saa > 1:
    seq=seq[:-(saa-1)] # removing the amino acids present before the amino acid from where the synthesis starts

dir = 'output/'
if not os.path.exists(dir):
    os.mkdir(dir)
os.chdir(dir) # changing current working directory to output folder
file = open(filename, 'w') # creating a new output file
file.write(datetime.datetime.now().strftime('%m-%d-%Y %I:%M:%S %p') + '\n')
file.write("The peptides sequence is " + seq + '\n')
file.close()
aan = len(seq)
ignore = ['*']
seq1 = Counter(x for x in seq if x not in ignore) # aa sorting
aan1 = len(seq1) # number of different amino acids and reagents
aa = [] # seq reversed for synthesis
a = [] # position
c = [] # coupling
d = [] # deprotection
ps.close()
    
aand = aan1 - (ports - 7)
if aand <= 0:
    paan = aan
    for n in range (1, paan+1):
        aa.append(seq[paan-n])
    positions(seq)
    presyn()
    syn()
else:
    if aand == 1 or aand == 2: 
        seqp1 = seq[5:aan+1]
        seqp2 = seq[0:5]
    if aand > 2 and aand < 10:
        seqp1 = seq[aan-ports+7:aan+1]
        seqp2 = seq[0:aan-ports+7]
    print("First part of the sequence to be synthesized is " + seqp1)
    print("Second part of the sequence to be synthesized is " + seqp2)
    file = open(filename, 'a')
    file.write("First part of the sequence to be synthesized is " + seqp1 + '\n')
    file.write("Second part of the sequence to be synthesized is " + seqp2 + '\n')
    file.close()
    print(" ")
    paan = len(seqp1)
    for n in range (1, paan+1):
        aa.append(seqp1[paan-n])
    positions(seqp1)
    print("First part of the sequence synthesis started")
    file = open(filename, 'a')
    file.write("First part of the sequence synthesis started" + '\n')
    file.close()
    presyn()
    syn()
    print("First part of the peptide synthesis done, amino acid/reagent lines will be cleaned")
    print("--------------------------------------------------------")
    print(" ")
    aalinecleaning()
    a = None
    c = None
    d = None
    aa = None
    a = []
    c = [] 
    d = [] 
    aa = []
    paan = len(seqp2)
    for n in range (1, paan+1):
        aa.append(seqp2[paan-n])
    positions(seqp2)
    print("Second part of the sequence synthesis started")
    file = open(filename, 'a')
    file.write("Second part of the sequence synthesis started" + '\n')
    file.close()
    syn()
    print("Second part of the peptide synthesis done")
    print("--------------------------------------------------------")
    print(" ")        
       
if fw == "y" or fw == "Y":
    finalwashing()
elif fw == "n" or fw == "N":
    print("Final washing skipped")
    print(" ")
else:
    fw1 = input("Input error in the sequence file. Do you want to perform final washing (y or n)? ")
    print(" ")
    if fw1 == "y" or fw1 == "Y":
        finalwashing()
    else:
        print("Final washing skipped")
        print(" ")
      
file = open(filename, 'a')
file.write("Peptide synthesis completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
file.close()

clean = input("Do you want to clean the amino acid/reagent lines (y or n)? ")
print(" ")
if clean == "y" or clean == "Y":
    print("--------------------------------------------------------")
    print("Amino acid/reagent lines cleaning")
    print("--------------------------------------------------------")
    print("Started at " + datetime.datetime.now().strftime('%I:%M:%S %p'))
    aalinecleaning()
    print("Completed at " + datetime.datetime.now().strftime('%I:%M:%S %p'))

print(" ")
print("PEPTIDE SYNTHESIS COMPLETED")
print("--------------------------------------------------------------------------------------------------")
# -------------------------------------------------------------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------------------------------------------------------------
