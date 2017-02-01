## PepSy
### An open-source automated peptide synthesizer based on Arduino and Python
Developed by Dr. Hariprasad Gali, Ph.D., Research Assistant Professor, Department of Pharmaceutical Sciences, College of Pharmacy, The University of Oklahoma Health Sciences Center, Oklahoma City, OK 73117.

Email address to report bugs: hgali@ouhsc.edu.
### Instructions
1. Install the following packages:
  1. Python (https://www.python.org/downloads/)
  2. pyserial (https://pypi.python.org/pypi/pyserial)
  3. pyFirmata (https://pypi.python.org/pypi/pyFirmata)
  4. Arduino Software (https://www.arduino.cc/en/Main/Software)
2. Upload standard firmata to the Arduino board from the File menu on the Arduino IDE, select Examples/Firmata/Standard Firmata and upload the file to the Arduino board.
3. Create folders named "sequence" and "output" within the same folder where PepSy.py and PepSy-manual.py scripts are saved.
4. Save device configuration file (config.txt) in the same folder where PepSy.py and PepSy-manual.py scripts are saved.
5. Create a sequence configuration file (see example templete.txt) for each run and save it in the "sequence" folder.
6. An output file is generated for each run and saved in the "output" folder.
7. PepSy.py script is written for operating the PepSy in a fully automatic mode.
8. PepSy-manual.py script is written for operating the PepSy in a fully manual mode and to clean amino acid/reagent lines.
9. Scripts were tested only with Python 3.5.0.
