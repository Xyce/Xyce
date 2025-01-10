%
% webwritenoproxy()
%
% This function mimics Matlab's webread() using Matlab's http class to 
% skip system proxy usage and force traffic to stay local
%

function response = webwritenoproxy( url, jsonMsg )
  % set up http options with proxy use off
  httpopt=matlab.net.http.HTTPOptions('UseProxy',0);
  method = matlab.net.http.RequestMethod.POST;
  mediaType = matlab.net.http.MediaType('application/json');
  contentType = matlab.net.http.field.ContentTypeField('text/plain')
  header = matlab.net.http.HeaderField('Media-Type', 'application/json', 'Content-Type', 'application/json') ;
  body = matlab.net.http.MessageBody(jsonMsg);
  
  %infos=containers.Map;
  request= matlab.net.http.RequestMessage(method, header, body );
  %request.field.ContentTypeField('text/plain');
  %request.Method = 'POST';
  %request.Body = jsonMsg;
  uri=matlab.net.URI(url);
  response=request.send(uri,httpopt);
end
