
This directory contains examples of calling Xyce from Python.

To use the examples here one must:

1. Install or build a copy of Xyce with shared library support on.  
   To see if an existing copy of Xyce has shared libraries look for
   libxycecinterface.dll under Windwos, libxycecinterface.so under Linux or
   libxycecinterface.dylib  under MacOS in the "bin" or "lib" directory
   where Xyce is installed (i.e. /usr/local/Xyce_7.8/lib/)  
   
2. Set the PYTHONPATH environment to where Xyce's "share" directory is located
   as in:  PYTHONPATH=/usr/local/Xyce_7.8/share/  This allows python to find 
   the python module xyce_interface.py  
   
3. Enter one of the example directories and run an examples. 
   cd runACircuit
   python runACircuit.py
   
   If needed, one can pass the location of the Xyce library directory to the 
   python program as
   
   python runACircuit.py /usr/local/Xyce_7.8/lib/  
   
   This can be necessary if xyce_interface.py hasn't had the library location 
   set properly during the installation process. 
   
   
   