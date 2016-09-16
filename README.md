## PepSy
### An open source automated peptide synthesizer based on Arduino and Python



Developed by Dr. Hariprasad Gali, Ph.D., Research Assistant Professor, Department of Pharmaceutical Sciences, College of Pharmacy, The University of Oklahoma Health Sciences Center, Oklahoma City, OK 73117.

Email address to report bugs: hgali@ouhsc.edu.

Tested only with Python 3.5.0.


Install the following packages:

1) Python (https://www.python.org/downloads/)

2) pyserial (https://pypi.python.org/pypi/pyserial)

3) pyFirmata (https://pypi.python.org/pypi/pyFirmata)

4) Arduino Software (https://www.arduino.cc/en/Main/Software)


Upload standard firmata to the Arduino board from the File menu on the Arduino IDE, select Examples/Firmata/Standard Firmata and upload the file to the Arduino board.

Create folders named "sequence" and "output" within the same folder where PepSy.py and PepSy-manual.py scripts are saved.

Device configuration file (config.txt) should be saved in the same folder where PepSy.py and PepSy-manual.py scripts are saved.

Sequence configuration file (see example templete.txt) should be saved in the "sequence" folder.

An output file is generated for each run and saved in the "output" folder.

PepSy.py script is for operating the PepSy in a fully automatic mode.

PepSy-manual.py script is for operating the PepSy in a fully manual mode and to clean amino acid/reagent line.
