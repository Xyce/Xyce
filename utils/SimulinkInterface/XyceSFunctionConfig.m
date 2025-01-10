
INIT FUNCTION

% Need to add the location of the simulink interface library to Matlabs search path.
%
% dumb way to do it for now.
addpath('/Users/rlschie/src/XyceDevelopment/BUILD/NormalCmake/utils/SimulinkInterface')
% 
% a better way would be it 
% 1. check for xyce_sfunction.mexmaci64 where the suffix will be platform dependent 
%    as in mexw64 (windows), mexa64 (linux) and mexmaci64 (MacOS)
% 2. if the file is not found then use pathtool to let the use select the path.



