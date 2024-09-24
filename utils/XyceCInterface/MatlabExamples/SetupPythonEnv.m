% have to tell the python environment within matlab where to find xyce_interface 
if count(py.sys.path,'/usr/local/Xyce_7.9/share') ==0 
  insert(py.sys.path,int32(0),'/usr/local/Xyce_7.9/share')
end


