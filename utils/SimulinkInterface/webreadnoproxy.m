%
% webreadnoproxy()
%
% This function mimics Matlab's webread() using Matlab's http class to 
% skip system proxy usage and force traffic to stay local
%

function response = webreadnoproxy( url )
  % set up http options with proxy use off
  httpopt=matlab.net.http.HTTPOptions('UseProxy',0);
  infos=containers.Map;
  request= matlab.net.http.RequestMessage;
  uri=matlab.net.URI(url);
  response=request.send(uri,httpopt);
end
