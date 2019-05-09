Xyce Parallel Electronic Simulator 

The binary distributions of Xyce reference a number of shared libraries that
are distributed along with the Xyce executable.  

The general command line usage to invoke Xyce is

    Xyce [Xyce options] netlistfile

Parallel builds of Xyce must be run through the mpirun or mpiexec program
provided by the OpenMPI installation:

    mpirun [mpirun options] Xyce [Xyce options] netlistfile

Both of the commands shown above assume that you have added the
directory in which Xyce's binaries were installed into your PATH
variable.