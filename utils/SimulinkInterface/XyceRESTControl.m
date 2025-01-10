
function XyceRESTControl

  % Create figure window
  mainFig = uifigure('Name', 'XyceRESTControl', 'CloseRequestFcn', @quitXyceRestControl, 'Position', [50, 100, 200, 300]);

  %layout 
  gl = uigridlayout(mainFig,[4 2]);
  gl.RowHeight = {'1x','1x','1x','1x'};
  gl.ColumnWidth = {'1x','5x'};

  % buttons for finding Xyce and Flask
  selXyceBt = uibutton(gl, 'Text', 'Locate Xyce', 'ToolTip', 'Locate Xyce on your system.', 'ButtonPushedFcn', @findXyce);
  selXyceBt.Layout.Row = 1;
  selXyceBt.Layout.Column = 2;
  selFlaskBt = uibutton(gl, 'Text', 'Locate Flask', 'ToolTip', 'Locate Python Flask on your system.', 'ButtonPushedFcn', @findFlask);
  selFlaskBt.Layout.Row = 2;
  selFlaskBt.Layout.Column = 2;

  % button to start/stop the running flask & Xyce process
  startStopBt = uibutton(gl, 'Text', 'Start', 'ToolTip', 'Start XyceRest Interface', 'Enable', 'off','ButtonPushedFcn', @toggleRunning);
  startStopBt.Layout.Row = 3;
  startStopBt.Layout.Column = 2;

  % button to quit
  quitBt = uibutton(gl, 'Text', 'Quit', 'ToolTip', 'Quit XyceRest Interface', 'ButtonPushedFcn', @quitXyceRestControl);
  quitBt.Layout.Row = 4;
  quitBt.Layout.Column = 2;
  
  % status lamps
  xyceStatusLamp = uilamp(gl, 'ToolTip', 'Indicates if Xyce is found and ready.');
  xyceStatusLamp.Layout.Row = 1;
  xyceStatusLamp.Layout.Column = 1;
  xyceStatusLamp.Color = 'red';

  flaskStatusLamp = uilamp(gl, 'ToolTip', 'Indicates if Flask is found and ready.');
  flaskStatusLamp.Layout.Row = 2;
  flaskStatusLamp.Layout.Column = 1;
  flaskStatusLamp.Color = 'red';
  
  runningStatusLamp = uilamp(gl, 'ToolTip', 'Indicates if Xyce REST interface is running.');
  runningStatusLamp.Layout.Row = 3;
  runningStatusLamp.Layout.Column = 1;
  runningStatusLamp.Color = 'red';
  
  

  % Store needed data in a structure in UserData.  This allows each UI element to access
  % the data as needed.  
  % The UI element StartStopButton is stored so it can be enabled if Xyce and Flask are found
  %
  mainFig.UserData = struct( ...
    'XyceLocation', '', ...
    'XyceRestFileLocation', '', ...
    'XyceRestFileLocationQuoted', '', ...
    'FlaskLocation', '', ...
    'XyceReady', false, ...
    'FlaskReady', false, ...
    'Running', false, ...
    'FlaskProcess', '', ...
    'StartStopBt', startStopBt, ...
    'XyceLamp', xyceStatusLamp, ...
    'FlaskLamp', flaskStatusLamp, ...
    'RunningLamp', runningStatusLamp );
  
  % initialize
  initializeState(mainFig)

end

%
function setExpectedRESTFileLoc(src)
  fig = ancestor(src,"figure","toplevel");
  % use fig.UserData.XyceLocation to get REST file location
  fig.UserData.XyceRestFileLocation = fig.UserData.XyceLocation;
  if (ispc)
    %display(fig.UserData.XyceLocation);
    fig.UserData.XyceRestFileLocation = replace(fig.UserData.XyceLocation, "bin\Xyce.exe", "share\XyceRest.py");
  elseif (ismac)
    fig.UserData.XyceRestFileLocation = replace(fig.UserData.XyceLocation, "bin/Xyce", "share/XyceRest.py");
  elseif (isunix)
    fig.UserData.XyceRestFileLocation = replace(fig.UserData.XyceLocation, "bin/Xyce", "share/XyceRest.py");
  end
  fig.UserData.XyceRestFileLocationQuoted = join( ['"',fig.UserData.XyceRestFileLocation, '"'], '');
end

% 
% function to let user find Xyce -- attached to Button
function findXyce(src,event)
  fig = ancestor(src,"figure","toplevel");
  [FileName, Path ] = uigetfile('*');
  fig.UserData.XyceLocation = strip(append(Path, FileName));
 
  setExpectedRESTFileLoc(src);

  % user selected a file, but we should check that we can 
  % find Xyce and XyceRest.py before assuming this is true.
  pyStatus = exist(fig.UserData.XyceRestFileLocation, 'file');
  if( (exist(fig.UserData.XyceLocation, 'file') == 2) && (exist(fig.UserData.XyceRestFileLocation, 'file') == 2))
    fig.UserData.XyceReady = true;
    fig.UserData.XyceLamp.Color = 'green';
  elseif( (exist(fig.UserData.XyceLocation, 'file') == 2) && (exist(fig.UserData.XyceRestFileLocation, 'file') == 0))
    fig.UserData.XyceReady = false;
    fig.UserData.XyceLamp.Color = 'yellow';
  else
    fig.UserData.XyceReady = false;
    fig.UserData.XyceLamp.Color = 'red';
  end

  if( fig.UserData.XyceReady && fig.UserData.FlaskReady )
    aUIButton = fig.UserData.StartStopBt;
    aUIButton.Enable = 'on';
  end
end

%
% function to let user find Flask -- attached to Button
function findFlask(src,event)
  fig = ancestor(src,"figure","toplevel");
  [FileName, Path ] = uigetfile('*');
  fig.UserData.FlaskLocation = append(Path, FileName);
  % user picked a file but we should do some sort of check to verify it is flask
  fig.UserData.FlaskReady = true;
  flaskLamp = fig.UserData.FlaskLamp;
  if( fig.UserData.FlaskReady )
    flaskLamp.Color = 'green';
  else
    flaskLamp.Color = 'red';
  end
  if( fig.UserData.XyceReady && fig.UserData.FlaskReady )
    aUIButton = fig.UserData.StartStopBt;
    aUIButton.Enable = 'on';
  end
end

%
% function to start and stop the XyceRest interface -- attached to Button
function toggleRunning(src,event)
  fig = ancestor(src,"figure","toplevel");
  if( ~fig.UserData.Running )
    result = startXyceRest( fig );
    if( result )
      fig.UserData.Running = true;
      aUIButton = fig.UserData.StartStopBt;
      aUIButton.Text = 'Stop';
      runningLamp = fig.UserData.RunningLamp;
      runningLamp.Color = 'green';
    end
  else
    result = stopXyceRest(fig );
    if( result )
      fig.UserData.Running = false;
      aUIButton = fig.UserData.StartStopBt;
      aUIButton.Text = 'Start';
      runningLamp = fig.UserData.RunningLamp;
      runningLamp.Color = 'red';
    end
  end
  
end


% function to check for a command on the users path
function valFound = findCommandOnUserPath( neededcmd )
  valFound = '';
  if (ispc)
    testcmd = append( 'where.exe ', neededcmd );
    [status, cmdout] = system(testcmd);
    if( status == 0 )
      valFound = strip(cmdout);
    end
  elseif (isunix || ismac)
    testcmd = append( 'which ', neededcmd );
    [status, cmdout] = system(testcmd);
    if( status == 0 )
      valFound = strip(cmdout);
    end
  end
end

% function to look for flask in OS dependent places.
function flaskFile = findFlaskOSSpecificLoc( )
  valFound = '';
  % on a mac look in ~/Library/Python/3.9/bin
  % under linux look in ~/.local/bin
  % windows \Users\Username\AppData\Local\Programs\Python\Python311 has python 
  % but packages are locally installed in a much more arcane path. 
  % "pip show flask" reports the exact package location so that might help.
  flaskFile='';
  if (ispc)
    %
    [status, cmdout] = system("pip show flask");
    locationStart = strfind(cmdout, 'Location:');
    locationEnd = strfind(cmdout, 'local-packages');
    location = cmdout(locationStart+10:locationEnd+13);
    pythonLongVer = regexp(location, 'Python\.\d\.\d\d', 'match');
    pythonLongVer = pythonLongVer{1};
    pythonVer = replace(pythonLongVer,'.', '');
    potentialFlaskLocation = [location, pythonVer, "Scripts", "flask.exe"];
    potentialFlaskLocation = join( potentialFlaskLocation, '\');
    %display('looking for flask in');
    %display(potentialFlaskLocation);
    ans = exist(potentialFlaskLocation);
    if( ans == 2)
      flaskFile = potentialFlaskLocation;
    end
  elseif (ismac)
    % look for local PIP install locaiton
    pipPyVer = strip(ls('~/Library/Python'));
    %[rows,cols] = size(pipPyVer);
    %pipPyVer = pipPyVer(1:(cols-1));
    usrName = getenv('USER');
    potentialFlaskLocation = ["/Users", usrName, "Library/Python", pipPyVer, "bin/flask"];
    potentialFlaskLocation = join(potentialFlaskLocation,'/');
    ans = exist(potentialFlaskLocation);
    if( ans == 2)
      flaskFile = potentialFlaskLocation;
    end
  elseif (isunix)
    % look for local PIP install locaiton
    usrName = getenv('USER');
    potentialFlaskLocation = ["/ascldap/users", usrName, ".local/bin/flask"];
    potentialFlaskLocation = join(potentialFlaskLocation,'/');
    ans = exist(potentialFlaskLocation);
    if( ans == 2)
      flaskFile = potentialFlaskLocation;
    end
  end
end

% functions to start and stop XyceRest interface
function val = startXyceRest( fig )
  val = true;
  flaskCommand = "";
  psCommand = "";

  % start flask and get the PID so we can stop it later.
  if (ispc)
    flaskCommand = join([fig.UserData.FlaskLocation, ' --app ', fig.UserData.XyceRestFileLocationQuoted, ' run &']);
    %('flask command');
    %display(flaskCommand)
    [ retcode, output ] = system(flaskCommand);
    psCommand = 'powershell -Command "Get-Process -Name flask | Format-List ID"';
    [ retcode, output ] = system(psCommand);
    lenOutput = length(output);
    fig.UserData.FlaskProcess = output(7:lenOutput-4);
    %display("flask process is:")
    %display(fig.UserData.FlaskProcess)
  elseif (ismac)
    flaskCommand = join([fig.UserData.FlaskLocation, " --app ", fig.UserData.XyceRestFileLocation, " run &"]);
    [ retcode, output ] = system(flaskCommand);
    psCommand = "ps -au $USER | grep flask | grep -v grep | cut -c7-11";
    [ retcode, output ] = system(psCommand);
    fig.UserData.FlaskProcess = output;
    %display("flask process is");
    %display(fig.UserData.FlaskProcess);
  elseif (isunix)
    % FLASK_APP=/fgs/rlschie/XyceDevelopment/INSTALL_SHARED/share/XyceRest.py ~/.local/bin/flask run
    flaskCommand = join(["FLASK_APP=", fig.UserData.XyceRestFileLocation, " ", fig.UserData.FlaskLocation, " run &"], '');
    [ retcode, output ] = system(flaskCommand);
    psCommand = "ps -au $USER | grep flask | grep -v grep | cut -c1-7";
    [ retcode, output ] = system(psCommand);
    fig.UserData.FlaskProcess = output;
    %display("flask process is");
    %display(fig.UserData.FlaskProcess);
  end
end

function val = stopXyceRest( fig )
  val = true;
  if (ispc)
    %killCommand = append('powershell.exe -Command "taskkill /F /IM flask.exe"');
    killCommand = append('powershell.exe -Command "Stop-Process -Force -Id ', fig.UserData.FlaskProcess, '"');
    %display("kill command is:");
    %display(killCommand);
    [retcode, output ] = system(killCommand);
  elseif( ismac)
    killCommand = append( 'kill -9 ', fig.UserData.FlaskProcess);
    [retcode, output ] = system(killCommand);
  elseif( isunix)
    killCommand = append( 'kill -9 ', fig.UserData.FlaskProcess);
    [retcode, output ] = system(killCommand);
  end
end


function initializeState( topFig )
  % check for Xyce
  xyceFound = findCommandOnUserPath( 'Xyce');
  if( strcmp(xyceFound, ''))
    display('Xyce not found.');
  else
    display('Xyce found');
    setExpectedRESTFileLoc(topFig);
    topFig.UserData.XyceLocation = xyceFound;
    topFig.UserData.XyceReady = true;
    if( (exist(topFig.UserData.XyceLocation, 'file') == 2) && (exist(topFig.UserData.XyceRestFileLocation, 'file') == 2))
      topFig.UserData.XyceReady = true;
      topFig.UserData.XyceLamp.Color = 'green';
      display('Xyce REST interface found')
    elseif( (exist(topFig.UserData.XyceLocation, 'file') == 2) && (exist(topFig.UserData.XyceRestFileLocation, 'file') == 0))
      topFig.UserData.XyceReady = false;
      topFig.UserData.XyceLamp.Color = 'yellow';
      display('Xyce REST interface not found.');
    else
      topFig.UserData.XyceReady = false;
      topFig.UserData.XyceLamp.Color = 'red';
    end
  end
  
  % check for Flask 
  flaskFound = findCommandOnUserPath( 'flask');
  if( strcmp(flaskFound, ''))
    % flask was not found on the path.  Look in OS dependent places.
    flaskFound = findFlaskOSSpecificLoc();
  end
  if( strcmp(flaskFound, ''))
    display('Flask not found.');
  else
    display('Flask found.');
    topFig.UserData.FlaskLocation = flaskFound;
    topFig.UserData.FlaskReady = true;
    topFig.UserData.FlaskLamp.Color = 'green';
  end
  
  if( topFig.UserData.XyceReady && topFig.UserData.FlaskReady )
    aUIButton = topFig.UserData.StartStopBt;
    aUIButton.Enable = 'on';
  end
  
end 
%
% function to close the main window
%
function quitXyceRestControl(src,event)
  fig = ancestor(src,"figure","toplevel");
  if(fig.UserData.Running)
    stopXyceRest(fig);
  end
  delete(fig)
end
