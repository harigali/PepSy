# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PepSy - An open-source peptide synthesizer
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Developed by Dr. Hariprasad Gali, Ph.D., Associate Professor of Research, Department of Pharmaceutical Sciences, College of Pharmacy, The University of Oklahoma Health Sciences Center, Oklahoma City, OK 73117.
# Email address to report bugs: hgali@ouhsc.edu.
# Tested only with Python 3.5.0
# Last update - July 23, 2020

# This script is written for synthesizing peptides using traditional fmoc chemistry. The synthesis conditions are optimized for 50 or 100 umol scale.
# This script includes ivDde deprotection, on-resin oxidation by Tl(CF3COO)3, and end capping with acetic anhydride.
# The synthesis can be paused before coupling of an amino acid.
# The synthesis conditions were only tested for the Rink Amide MBHA resin and with the Wang resin and Cl-Trt resin with the first amino acid already coupled.

# Create folders named 'sequence' and 'output' within the same folder where PepSy.py and PepSy-manual.py scripts are saved.
# Save device configuration file (config.txt) in the same folder where PepSy.py and PepSy-manual.py scripts are saved.
# Create a sequence configuration file (see example templete.txt) for each run and save it in the 'sequence' folder.
# An output file is generated for each run and stored in the 'output' folder.

# Only Arduino UNO digital pins are used.
# COM port numbers of VICI stream selector valve (ps) and Arduino UNO (board) needs to be updated in the device configuration file according to their current assignment on the PC.

# VICI CHEMINERT low pressure 24 stream selector valve
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
# Amino acid priming volumes are optimized for tube lengths (~11.6 microliter per inch for 0.03' ID tubing) - ~15 inch from aa to ps, ~11 inch from ps to pump, and ~20 inch pump to resin.
# Tubing volume variables are len1 (aa to ps),  len2 (ps to pump), and  len3 (pump to resin) depds on tubing id and lengths, which needs to be updated in the device configuration file when replacing the tubing with different length or id.
# Solution concentrations to be used: amino acids - 0.33M, HBTU - 0.33M, HOBT - 0.66M, DIPEA - 1.32M, piperidine - 20%, hydrazine - 2%, acetic anhydride/pyridine - 2.5M/2.5M, and Tl(CF3COO)2 - 0.05M.
# If the solenoid micro pump is replaced with a different pump internal volume (piv), it needs to be updated in the device configuration file.

# Uppercase alphabets are used for L amino acids.
# Lowercase alphabets are used for D amino acids.
# '<', '>', '-', '+', and '=' are used for N-methyl amino acids.
# '#' is used for a manually added amino acid. 
# '3', '4', '5', '6', and '8' are used for beta-alanine, 4-aminobutanoic acid, 5-aminovaleric acid, 6-aminohexanoic acid, and 8-aminooctanoic acid linker and place the solution in the position assigned to '3', '4', '5', '6', or '8'  respectively
# 'X' and 'B' are used for  H2N-PEG2-COOH and H2N-PEG3-COOH linker respectively
# 'J', '1', '2', '7', or '9' are used for a linker (other than beta-alanine, 4-aminobutanoic acid, 5-aminovaleric acid, 6-aminohexanoic acid, 8-aminooctanoic acid, H2N-PEG2-COOH, H2N-PEG3-COOH), an unusual amino acid or any molecule that requires both coupling and fmoc deprotection and place the solution in the position assigned to 'J', '1', '2', '7', or '9' respectively.
# 'Z' is used for tris t-butyl protected DOTA and place the DOTA solution in the position assigned to 'Z'.
# 'U' or 'O' are used for a chelator (other than tris t-butyl protected DOTA), an unusual amino acid, or any molecule that requires only coupling and place the solution in the position assigned to  'U' or 'O' respectively.
# '*' is used for pausing the synthesis.
# '!' is used for ivDde deprotection and place the hydrazine solution in the position assigned to '!'.
# '@' is used for onresin oxidation and add the thallium solution manually.
# '$' is used for endcapping and place the acetic anhydride solution in the position assigned to '$'.
# -------------------------------------------------------------------------------------------------------------------------------------------

# Imports
from os import path, mkdir, chdir
from time import sleep
from datetime import datetime
from configparser import ConfigParser
from pyfirmata import Arduino, util
from collections import Counter
import serial
# -------------------------------------------------------------------------------------------------------------------------------------------

# Functions
def timestamp():
    timestamp = datetime.now().strftime('%I:%M:%S %p')
    return timestamp

def filewrite(info):
    file = open(filename, 'a')
    file.write(info + '\n')
    print(info)
    file.close()
    
def positions(p):
    mwdict = {'A':329.36, 'C':585.72, 'D':411.45, 'E':425.48, 'F':387.44, 'G':297.31, 'H':619.72, 'I':353.42, 'K':468.2, 'L':353.42, 'M':371.45, 'N':596.68, 'P':337.38, 'Q':610.71, 'R':648.78, 'S':383.44,
              'T':379.48, 'V':339.39, 'W':526.59, 'Y':459.54, '3':311.3, '4':325.4, '5':339.4, '6':353.3, '8':381.5, 'X':385.42, 'B':429.47, 'Z':572.74} # molecular weight of standard fmoc-protected amino acids
    paap = [] # positions for different amino acids and reagents
    pseq = Counter(x for x in p if x not in ignore) # amino acids and reagents sorting
    paan1 = len(pseq) # number of different amino acids and reagents
    for n in range (2, paan+1):
        if aa[n-2] == 'P' or aa[n-2] in ('P', '<', '>', '+', '-', '='):
            pseq[aa[n-1]] += 1               
    paak = list(pseq.keys())
    paav = list(pseq.values())  
    if pa.upper() == 'Y':
        filewrite('Place amino acid/reagent solutions with required volumes in the positions shown below')
        print(' ')
        filewrite('---------------------------------------------------------------------------------------------------')
        filewrite('S. No.' + '\t' + 'Amino acid' + '\t' + 'Position' + '\t' + 'Solution volume' + '\t\t' + 'Amino acid weight' + '\t' + 'DMF volume')
        filewrite('---------------------------------------------------------------------------------------------------')
        for n in range (1, paan1+1):
            at = n % (ports - 7)
            if at == 0:
                pos = ports
            else:
                pos = 7 + at
            paap.append(pos)
            if paak[n-1] == '!':
                vol = ((len1+len2)/1000+ss*2)*paav[n-1]
            elif paak[n-1] == '@':
                vol = ((len1+len2)/1000+ss*4)*paav[n-1]
            elif paak[n-1] == '$':
                vol = ((len1+len2)/1000+ss*1)*paav[n-1]
            else:
                vol = ((len1+len2)/1000+ss*0.5)*paav[n-1]  
            try:
                mw = mwdict[paak[n-1].upper()]
            except:
                mw = 0
            wt = vol*mw*0.33 # conc of aa solution is 0.33M
            dmf = vol*1000-wt
            filewrite(str(n) + '\t' + paak[n-1] + '(' + str(paav[n-1]) + ')' + '\t\t' + str(paap[n-1]) + '\t\t' + str("{:.1f}".format(vol)) + ' ml' + '\t\t\t' + str("{:.0f}".format(wt)) + ' mg' + '\t\t\t' + str("{:.0f}".format(dmf)) + ' ul')
        filewrite('---------------------------------------------------------------------------------------------------')
        print(' ')
    else:
        for n in range (1, paan1+1):
            pos = synconfig.getint('Positions', paak[n-1])
            paap.append(pos)
    filewrite('-----------------------------------------------------------------------------')
    filewrite('S. No.' + '\t' + 'Amino acid' + '\t' + 'Position' + '\t' + 'Coupling' + '\t\t' + 'Deprotection')
    filewrite('-----------------------------------------------------------------------------')
    for n in range (1, paan+1):
        if aa[n-1] in ('*', '@', '#'):
            a.append(1)
        for m in range (1, paan1+1):
            if aa[n-1] == paak[m-1]:
                a.append(paap[m-1])
        if aa[n-1] == '!':
            c.append('ivdde')
        elif aa[n-1] == '@':
            c.append('oxidation')
        elif aa[n-1] == '$':
            c.append('encapping')
        elif aa[n-1] == '*':
            c.append('pause')
        elif aa[n-1] == '#':
            c.append('manual')
        else:
            c.append('single') # default coupling
        if n > 1:
            if aa[n-2] in ('P', '<', '>', '+', '-', '='):
                c.pop() # deletes the last element from the list
                c.append('double') # double coupling if previous aa is P or any aa represented by <, >, +, -, or =
        if aa[n-1] in ('*', '!', '@', '$', 'Z', 'U', 'O'):
            d.append('none')
        else:
            d.append('fmoc') # default deprotection
        if saa > 1:
            n1 = n + saa - 1
        else:
            n1 = n
        filewrite(str(n1) + '\t' + aa[n-1] + '\t\t' + str(a[n-1]) + '\t\t' + c[n-1] + '\t\t\t' + d[n-1])
    filewrite('-----------------------------------------------------------------------------')
    dmfn = 0 # number of dmf washings
    cn = 0 # number of couplings
    dn = 0 # number of fmoc deprotections
    for n in range (1, paan+1):
        if aa[n-1] != '*':
            dmfn += 10 # 5 for coupling and 5 for deprotection
        if aa[n-1] in ('@', '!', '$', 'Z', 'U', 'O'):
            dmfn -= 5 # 5 for deprotection removed
        if aa[n-1] not in ('*', '@', '!', '$'):
            cn += 1
        if aa[n-1] in ('p', 'P', '<', '>', '-', '+', '='):
            cn += 2 # double coupling
        if aa[n-1] not in ('*', '@', '!', '$', 'Z', 'U', 'O'):
            dn += 1
    dmfvol = 5+11.5+(dmfn*2) # 5 ml for extra, 0.5 ml for initial priming/swelling/fmocdeprotection and 2 ml for each washing
    dcmvol = 10+11.5 # 10 ml for extra, 11.5 ml for initial priming, swelling, and final washing    
    hbtvol = 2+0.5+((len2/1000)+0.5)*cn # 2 ml for extra, 0.5 ml for initial priming, and (len2 + 0.5) for each coupling    
    hobvol = 2+0.5+((len2/1000)+0.26)*cn # 2 ml for extra, 0.5 ml for initial priming, and (len2 + 0.26) for each coupling    
    dipvol = 2+0.5+((len2/1000)+0.26)*cn # 2 ml for extra, 0.5 ml for initial priming, and (len2 + 0.26) for each coupling    
    pipvol = 2+0.5+((len2+len3)/1000+2)*dn # 2 ml for extra, 0.5 ml for initial priming, and (len2 + len3 + 2) for each fmoc deprotection    
    print(' ')
    filewrite('Volume of DMF = ' + str("{:.1f}".format(dmfvol)) + ' ml')
    filewrite('Volume of DCM = ' + str("{:.1f}".format(dcmvol)) + ' ml')
    filewrite('Volume of HBTU solution = ' + str("{:.1f}".format(hbtvol)) + ' ml')
    filewrite('Volume of HOBt solution = ' + str("{:.1f}".format(hobvol)) + ' ml')
    filewrite('Volume of DIPEA solution = ' + str("{:.1f}".format(dipvol)) + ' ml')
    filewrite('Volume of Piperidine solution  = ' + str("{:.1f}".format(pipvol)) + ' ml')
    print(' ')
    print('Check the positions, couplings, and deprotections are correct')
    print(' ')
    print('Check the nitrogen gas pressure, if it is not ~2 psi then adjust the pressure.')
    print(' ')
    print('Check the levels of all the reagents, if any of them is not enough then add.')
    print(' ')
    input('If you are ready, press ENTER to continue')
    
def pspos(p): # p is stream selector position (integer)
    if p == 1:
        position = 'HM\r'
    else:
        position = 'GO%d\r' % (p)
    ps.open()
    ps.write(position.encode())
    ps.close()
    
def pumpon(v): # v is volume (integer) to be pumped in microliters
    for p in range(0, v, piv): 
        pump.write(1)
        sleep(0.25)
        pump.write(0)
        sleep(0.25)

def presyn():
    initialization()
    if pr.upper() == 'Y':
        priming()
    elif pr.upper() == 'N':
        filewrite('Priming skipped')
        print(' ')
    else:
        pr1 = input('Input error in the sequence file. Do you want to perform priming (y or n)? ')
        print(' ')
        if pr1.upper() == 'Y':
            priming()
        else:
            filewrite('Priming skipped')
            print(' ')
    if sw.upper() == 'Y':
        swelling()
    elif sw.upper() == 'N':
        filewrite('Swelling skipped')
        print(' ')
    else:
        sw1 = input('Input error in the sequence file. Do you want to perform swelling (y or n)? ')
        print(' ')
        if sw1.upper() == 'Y':
            swelling()
        else:
            filewrite('Swelling skipped')
            print(' ')
    filewrite('Peptide synthesis started at ' + timestamp())
    print(' ')
    if dp.upper() == 'Y':
        fmocdeprotection()
    elif dp.upper() == 'N':
        print('Initial fmoc deprotection skipped')
        print(' ')
    else:
        dp1 = input('Input error in the sequence file. Do you want to perform initial fmoc deprotection (y or n)? ')
        print(' ')
        if dp1.upper() == 'Y':
            fmocdeprotection()
        else:
            filewrite('Initial fmoc deprotection skipped')
            print(' ')
        
def syn():
    for n in range(1, paan+1):
        if saa > 1:
            n1 = n + saa - 1
        else:
            n1 = n
        if aa[n-1] == '*':
            filewrite('Pause')
        elif aa[n-1] == '@':
            filewrite('On-resin oxidation')
        elif aa[n-1] == '!':
            filewrite('ivDde deprotection')
        elif aa[n-1] == '$':
            filewrite('Endcapping')
        else:
            filewrite('Amino acid: ' + str(n1) + ' (' + aa[n-1] + ')')
        if c[n-1] == 'single' or c[n-1] == 'manual': 
            coupling(n-1)
        elif c[n-1] == 'double':
            doublecoupling(n-1)
        elif c[n-1] == 'oxidation':
            onresinoxidation()
        elif c[n-1] == 'endcapping':
            endcapping(n-1) 
        elif c[n-1] == 'ivdde':
            ivddedeprotection(n-1) 
        elif c[n-1] == 'pause':
            pause()
        else:
            filewrite('No coupling')
            print(' ')
        if d[n-1] == 'fmoc':
            fmocdeprotection()

def washing():
    reagent.write(1)
    pspos(2)
    pumpon(2000) # addition of 2 ml DMF
    pspos(1)
    sleep(1)
    reagent.write(0)
    n2.write(1)
    waste.write(1)
    vent.write(1)
    sleep(60) # draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
     
def initialization():
    filewrite('Initialization started at ' + timestamp())
    pspos(1)
    n2.write(0)
    vent.write(0)
    reagent.write(0)
    waste.write(0)
    prime.write(0)
    pump.write(0)
    filewrite('Completed at ' + timestamp())
    print(' ')

def priming():
    filewrite('Priming started at ' + timestamp())
    prime.write(1)
    for p in range(4,8):
        pspos(p)
        pumpon(len1+len2)
        pspos(1)
        sleep(1)
    for p in range(2,4):
        pspos(p)
        pumpon(1000)
        pspos(1)
        sleep(1)
    prime.write(0)
    filewrite('Completed at ' + timestamp())
    print(' ')

def swelling():
    print('Swelling started at ' + timestamp())
    print('Adding solvents')
    reagent.write(1)
    pspos(3)
    pumpon(1000-len3) # addition of 1 ml DCM, 768 from pumping, 232 from tubing
    pspos(1)
    sleep(1)
    pspos(2)
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(1000+len3) # addition of 1 ml DMF, 1000 DMF, 232 residual DCM
    pspos(1)
    reagent.write(0)
    n2.write(1)
    print('15 min swelling')
    sleep(900) # 15 min swelling
    print('Draining solvents')
    waste.write(1)
    vent.write(1)
    sleep(30) # draining
    n2.write(0)
    waste.write(0)
    vent.write(0)
    print('Washing')
    washing() # washing
    filewrite('Completed at ' + timestamp())
    print(' ')
    
def coupling(n): # S. No. (integer) of the amino acid
    filewrite('Coupling (single) started at ' + timestamp())
    if a[n] == 1:
        input('Synthesis paused, add amino acid solution to the reactor manually, and press ENTER to continue')
    else:
        aapos = a[n]
        pspos(aapos)
        prime.write(1)
        filewrite('Amino acid position on PS is ' + str(aapos))
        print('Priming amino acid ' + aa[n])
        pumpon(len1+len2) # amino acid line priming - aa to ps to pump
        prime.write(0)
        reagent.write(1)
        pumpon(len3) # amino acid line priming - pump to resin
        reagent.write(0)
        waste.write(1)
        vent.write(1)
        n2.write(1)
        sleep(10)
        n2.write(0)
        vent.write(0)
        waste.write(0)
        reagent.write(1)
        pumpon(ss*500-len3) # addition of 0.5 ml amino acid solution
        pspos(1)
        sleep(1)
    print('Adding reagents')
    pspos(5)
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(ss*260) # addition of 0.26 ml DIPEA solution
    pspos(1)
    sleep(1)
    pspos(6)
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(ss*260) # addition of 0.26 ml HOBT solution
    pspos(1)
    sleep(1)
    pspos(7)
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(ss*500) # addition of 0.5 ml HBTU solution
    pspos(1)
    reagent.write(0)
    pspos(2)
    prime.write(1)
    reagent.write(0)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing pump to resin
    pspos(1)
    reagent.write(0)
    n2.write(1)
    print('60 min coupling')
    sleep(3600) # 60 min coupling
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    n2.write(0)
    waste.write(0)
    vent.write(0)
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        washing()
    filewrite('Completed at ' + timestamp())
    print(' ')

def doublecoupling(n): # S. No. (integer) of the amino acid
    filewrite('Coupling (double) started at ' + timestamp())
    for d in range(1,3):
        aapos = a[n]
        pspos(aapos)
        prime.write(1)
        if d == 1:
            filewrite('Amino acid position on PS is ' + str(aapos))
            print('Priming amino acid ' + aa[n])
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
        sleep(10)
        n2.write(0)
        vent.write(0)
        waste.write(0)
        print('Adding reagents')
        reagent.write(1)
        pumpon(ss*500-len3) # addition of 0.5 ml amino acid solution
        pspos(1)
        sleep(1)
        pspos(5)
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # removing previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(ss*260) # addition of 0.26 ml DIPEA solution
        pspos(1)
        sleep(1)
        pspos(6)
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # removing previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(ss*260) # addition of 0.26 ml HOBT solution
        pspos(1)
        sleep(1)
        pspos(7)
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # removing previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(ss*500) # addition of 0.5 ml HBTU solution
        pspos(1)
        reagent.write(0)
        pspos(2)
        prime.write(1)
        reagent.write(0)
        pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
        prime.write(0)
        reagent.write(1)
        pumpon(len3) # DMF to add previous reagent leftover in the tubing
        reagent.write(0)
        pspos(1)
        n2.write(1)
        print('60 min coupling')
        sleep(3600) # 60 min coupling
        waste.write(1)
        vent.write(1)
        print('Draining reagents')
        sleep(30)# draining
        n2.write(0)
        waste.write(0)
        vent.write(0)
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        washing()
    filewrite('Completed at ' + timestamp())
    print(' ')

def fmocdeprotection():
    filewrite('fmoc deprotection started at ' + timestamp())
    print('Adding reagents')
    pspos(4)
    sleep(1)
    prime.write(1)
    pumpon(len2) # removing previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*1000) # addition of 1 ml piperidine solution
    pspos(1)
    n2.write(1)
    reagent.write(0)
    print('10 min first round deprotection')
    sleep(600) # 10 min first round deprotection
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    print('Adding reagents')
    reagent.write(1)
    pspos(4)
    pumpon(ss*1000) # addition of 1 ml piperidine solution
    pspos(1)
    n2.write(1)
    reagent.write(0)
    print('20 min second round deprotection')
    sleep(1200) # 20 min second round deprotection
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    pspos(2)
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing
    pspos(1)
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        washing()
    filewrite('Completed at ' + timestamp())
    print(' ')
    
def finalwashing():
    filewrite('Final washing tarted at ' + timestamp())
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        reagent.write(1)
        pspos(3)
        pumpon(2000) # addition of 2 ml DCM
        pspos(1)
        sleep(1)
        reagent.write(0)
        n2.write(1)
        waste.write(1)
        vent.write(1)
        sleep(60) # draining
        vent.write(0)
        waste.write(0)
        n2.write(0)
    filewrite('Completed at ' + timestamp())
    print(' ')
    filewrite('Drying started at ' + timestamp())
    drying()
    filewrite('Completed at ' + timestamp())
    print(' ')

def pause():
    reagent.write(1)
    pspos(2)
    pumpon(1000) # addition of 1 ml DMF
    pspos(1)
    reagent.write(0)
    input('Synthesis paused, Press ENTER to continue')
    print(' ')
    n2.write(1)
    waste.write(1)
    vent.write(1)
    sleep(15) # draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    
def aalinecleaning():
    m = int(input('Enter the starting position on ps '))
    n = int(input('Enter the ending position on ps '))
    print(' ')
    input('Insert all amino acid/reagent lines in DMF and then press ENTER to continue')
    print(' ')
    for o in range(m,n+1):
        print('Cleaning line ', o)
        pspos(o)
        prime.write(1)
        pumpon(500+len1+len2)
        prime.write(0)
    pspos(1)
    print(' ')
    print('Remove amino acid/reagent lines from DMF and clean the exterior with acetone or isopropyl alcohol wipe')
    print(' ')
    print('Amino acid/reagent lines cleaning completed')
    print(' ')
    
def drying():
    n2.write(1)
    vent.write(1)
    waste.write(1)
    sleep(1800) # 30 min drying
    n2.write(0)
    vent.write(0)
    waste.write(0)

def ivddedeprotection(n): # S. No. (integer) of the reagent
    filewrite('ivDde deprotection started at ' + timestamp())
    print('Adding reagents')
    aapos = a[n]
    pspos(aapos)
    filewrite('ivDde position on PS is ' + str(aapos))
    sleep(1)
    prime.write(1)
    pumpon(len1+len2) # removing previous reagent from tubing between aa to ps to pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*1000) # addition of 1 ml hydrazine solution
    pspos(1)
    n2.write(1)
    reagent.write(0)
    print('10 min first round deprotection')
    sleep(600) # 10 min first round deprotection
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    print('Adding reagents')
    reagent.write(1)
    pspos(aapos)
    pumpon(ss*1000) # addition of 1 ml hydrazine solution
    pspos(1)
    n2.write(1)
    reagent.write(0)
    print('20 min second round deprotection')
    sleep(1200) # 20 min second round deprotection
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    pspos(2)
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing
    pspos(1)
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        washing()
    filewrite('Completed at ' + timestamp())
    print(' ')

def onresinoxidation():
    filewrite('Onresin oxidation started at ' + timestamp())
    input('Synthesis paused, add Tl(CF3COO)3 solution to the reactor manually, and press ENTER to continue')
    n2.write(1)
    print('60 min first round oxidation')
    sleep(3600) # 60 min first round oxidation
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    input('Synthesis paused, add Tl(CF3COO)3 solution to the reactor manually, and press ENTER to continue')
    n2.write(1)
    print('60 min second round oxidation')
    sleep(3600) # 60 min second round oxidation
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    vent.write(0)
    waste.write(0)
    n2.write(0)
    pspos(2)
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len2) # DMF to add previous reagent leftover in the tubing
    pspos(1)
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        washing()
    filewrite('Completed at ' + timestamp())
    print(' ')

def endcapping(n): # S. No. (integer) of the reagent
    filewrite('Endcapping started at ' + timestamp())
    print('Adding reagents')
    aapos = a[n]
    pspos(aapos)
    filewrite('Acetic anhydride position on PS is ' + str(aapos))
    sleep(1)
    prime.write(1)
    pumpon(len1+len2) # removing previous reagent from tubing between aa to ps to pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # removing DMF leftover in the tubing
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    reagent.write(1)
    pumpon(ss*1000) # addition of 1 ml acetic anhydride solution
    pspos(1)
    n2.write(1)
    reagent.write(0)
    print('30 min end capping')
    sleep(1800) # 30 min end capping
    waste.write(1)
    vent.write(1)
    print('Draining reagents')
    sleep(30)# draining
    waste.write(0)
    vent.write(0)
    n2.write(0)
    pspos(2)
    prime.write(1)
    pumpon(len2) # DMF to remove previous reagent from tubing between ps and pump
    prime.write(0)
    reagent.write(1)
    pumpon(len3) # DMF to add previous reagent leftover in the tubing
    pspos(1)
    reagent.write(0)
    waste.write(1)
    vent.write(1)
    n2.write(1)
    sleep(10)
    n2.write(0)
    vent.write(0)
    waste.write(0)
    for w in range(1,6): # 5 times washing
        print('Washing ' + str(w))
        washing()
    filewrite('Completed at ' + timestamp())
    print(' ')
# -------------------------------------------------------------------------------------------------------------------------------------------

# Main
devconfig = ConfigParser()
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

print(' ')
print('--------------------------------------------------------------------------------------------------')
print('                                           PepSy                                                  ')
print('--------------------------------------------------------------------------------------------------')
print(' ')

seqfile = input('Enter the sequence configuration file name ')
print(' ')
if seqfile == '':
    q1 = input('Do you want to terminate the run (y or n)? ')
    print(' ')
    if q1.upper() == 'Y':
        ps.close()
        exit()
    else:
        seqfile = input('Enter the sequence configuration file name ')
        print(' ')
filename = seqfile + datetime.now().strftime('-%Y-%m-%d-') +  'out.txt'
seqfile = 'sequence/' + seqfile + '.txt'

synconfig = ConfigParser()
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
if not path.exists(dir):
    mkdir(dir)
chdir(dir) # changing current working directory to output folder
file = open(filename, 'w') # creating a new output file
file.close()
filewrite(datetime.now().strftime('%m-%d-%Y %I:%M:%S %p'))
filewrite('The peptides sequence not including any amino acid already present on the resin is ' + seq + '\n')
print(' ')
aan = len(seq)
ignore = ['*', '@', '#']
seq1 = Counter(x for x in seq if x not in ignore) # aa sorting
aan1 = len(seq1) # number of different amino acids and reagents
aa = [] # seq reversed for synthesis
a = [] # position
c = [] # coupling
d = [] # deprotection
ps.close()
    
if aan1 <= ports - 7: # checking whether reqired number of ports are available to accommodate all the amino acids and reagents, if not the peptide sequence will be split in to two parts
    paan = aan
    for n in range (1, paan+1):
        aa.append(seq[paan-n])
    positions(seq)
    presyn()
    syn()
else:
    i = 1
    while True:
        seqtemp = seq[i:aan+1]
        seqtemp1 = Counter(x for x in seqtemp if x not in ignore) # aa sorting
        aantemp = len(seqtemp1) # number of different amino acids and reagents    
        if aantemp <= ports - 7:
            seqp1 = seq[i:aan+1]
            seqp2 = seq[0:i]
            break
        i = i + 1
    filewrite('First part of the sequence to be synthesized is ' + seqp1)
    filewrite('Second part of the sequence to be synthesized is ' + seqp2)
    print(' ')
    paan = len(seqp1)
    for n in range (1, paan+1):
        aa.append(seqp1[paan-n])
    positions(seqp1)
    filewrite('First part of the sequence synthesis started')
    presyn()
    syn()
    filewrite('First part of the peptide synthesis done, amino acid/reagent lines will be cleaned')
    print(' ')
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
    filewrite('Second part of the sequence synthesis started')
    syn()
    filewrite('Second part of the peptide synthesis done')
    print(' ')        
       
if fw.upper() == 'Y':
    finalwashing()
elif fw.upper() == 'N':
    filewrite('Final washing skipped')
    print(' ')
else:
    fw1 = input('Input error in the sequence file. Do you want to perform final washing (y or n)? ')
    print(' ')
    if fw1.upper() == 'Y':
        finalwashing()
    else:
        filewrite('Final washing skipped')
        print(' ')
      
clean = input('Do you want to clean the amino acid/reagent lines (y or n)? ')
print(' ')
if clean.upper() == 'Y':
    filewrite('Amino acid/reagent lines cleaning started at ' + timestamp())
    aalinecleaning()
    filewrite('Completed at ' + timestamp())

print(' ')
filewrite('Peptide synthesis completed at ' + timestamp())
# -------------------------------------------------------------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------------------------------------------------------------
