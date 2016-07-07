# PepSy
### An open source peptide synthesizer


Utilizes Arduino Uno and python

Developed by Dr. Hariprasad Gali, Ph.D., Research Assistant Professor, Department of Pharmaceutical Sciences, The University of Oklahoma Health Sciences Center.

Email address to report bugs is hgali@ouhsc.edu.

Tested only with Python 3.5.0

Install pyFirmata.

Upload standard firmata to the Arduino Uno from the File menu on the Arduino IDE, select Examples/Firmata/Standard Firmata and upload the file to the Arduino.

Create folders named "sequence" and "output" within the same folder where PepSy.py and PepSy-manual.py scripts are saved.

Device configuration file (config.txt) should be saved in the same folder as this script.

Sequence configuration file (see example templete.txt) should be saved in the "sequence" folder.

An output file is generated for each run and saved in the "output" folder.
