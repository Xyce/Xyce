This directory contains examples of calling Xyce from Matlab/Simulink.

The Matlab/Simulink interface uses Xyce's Python interface from Matlab.  Thus,
one should ensure that they are using a version of Python compatible with
Matlab (see https://www.mathworks.com/support/requirements/python-compatibility.html).

Additionally, one should test the Python examples to ensure the Xyce - Python 
layer is installed and working correctly.  See the directory above this one 
labeled PythonExamples.

Finally the Matlab script in the Matlab Examples directory will need to edited 
for the location of your local installation of Xyce's share directory where 
the python file xyce_interface.py is located.

See the document Application Note: Using Xyce From Python, Matlab and Simulink for
more details.