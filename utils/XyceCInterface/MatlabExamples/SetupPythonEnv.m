% have to tell the python environment within matlab where to find xyce_interface 
if count(py.sys.path,'/Users/rlschie/src/XyceDevelopment/INSTALL/share') ==0 
  insert(py.sys.path,int32(0),'/Users/rlschie/src/XyceDevelopment/INSTALL/share')
end


