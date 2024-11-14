% have to tell the python environment within matlab where to find xyce_interface 
if count(py.sys.path,'@CPACK_PACKAGING_INSTALL_PREFIX@/share') ==0 
  insert(py.sys.path,int32(0),'@CPACK_PACKAGING_INSTALL_PREFIX@/share')
end
if count(py.sys.path,'@CMAKE_INSTALL_PREFIX@/share') ==0 
  insert(py.sys.path,int32(0),'@CMAKE_INSTALL_PREFIX@/share')
end

