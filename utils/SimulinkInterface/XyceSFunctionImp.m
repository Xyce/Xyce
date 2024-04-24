function XyceSFunctionImp(block)
%XyceSFunctionImp A Xyce mixed signal connection  for a MATLAB S-Function
%   A matlab s-function which connects Simulink to a Xyce instance 
%   through Xyces mixed signal interface. 
  
%
% The setup method is used to setup the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.  
%   
setup(block);
  
%endfunction

% Function: setup ===================================================
% Abstract:
%   Set up the S-function block's basic characteristics such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
% 
%   Required         : Yes
%   C MEX counterpart: mdlInitializeSizes
%
function setup(block)

  % Register the number of ports.
  block.NumInputPorts  = block.DialogPrm(3).Data;
  block.NumOutputPorts = block.DialogPrm(4).Data;
  
  % Set up the port properties to be inherited or dynamic.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  for idx = 1:1:block.NumInputPorts
    % Override the input port properties.
    block.InputPort(idx).DatatypeID  = 0;  % double
    block.InputPort(idx).Complexity  = 'Real';
  end
  
  for idx = 1:1:block.NumOutputPorts
    % Override the output port properties.
    block.OutputPort(idx).DatatypeID  = 0; % double
    block.OutputPort(idx).Complexity  = 'Real';
  end

  % Register the parameters.
  block.NumDialogPrms     = 6;
  block.DialogPrmsTunable = {'Nontunable','Nontunable','Nontunable','Nontunable','Nontunable','Nontunable'};
  
  %inputPortNames = block.DialogPrm(5).Data;
  %outputPortNames = block.DialogPrm(6).Data;
  %display('In setup inputPorts Data');
  %display(block.DialogPrm(5).Data);
  %display('In setup outputPorts Data');
  %display(block.DialogPrm(6).Data);
  
  % Set up the continuous states.
  block.NumContStates = 1;

  % Register the sample times.
  %  [0 offset]            : Continuous sample time
  %  [positive_num offset] : Discrete sample time
  %
  %  [-1, 0]               : Inherited sample time
  %  [-2, 0]               : Variable sample time
  block.SampleTimes = [0 0];
  
  % -----------------------------------------------------------------
  % Options
  % -----------------------------------------------------------------
  % Specify if Accelerator should use TLC or call back to the 
  % MATLAB file
  block.SetAccelRunOnTLC(false);
  
  % Specify the block's operating point compliance. The block operating 
  % point is used during the containing model's operating point save/restore)
  % The allowed values are:
  %   'Default' : Same the block's operating point as of a built-in block
  %   'UseEmpty': No data to save/restore in the block operating point
  %   'Custom'  : Has custom methods for operating point save/restore
  %                 (see GetOperatingPoint/SetOperatingPoint below)
  %   'Disallow': Error out when saving or restoring the block operating point.
  block.OperatingPointCompliance = 'Default';
  
  % -----------------------------------------------------------------
  % The MATLAB S-function uses an internal registry for all
  % block methods. You should register all relevant methods
  % (optional and required) as illustrated below. You may choose
  % any suitable name for the methods and implement these methods
  % as local functions within the same file.
  % -----------------------------------------------------------------
   
  % -----------------------------------------------------------------
  % Register the methods called during update diagram/compilation.
  % -----------------------------------------------------------------
  
  % 
  % CheckParameters:
  %   Functionality    : Called in order to allow validation of the
  %                      block dialog parameters. You are 
  %                      responsible for calling this method
  %                      explicitly at the start of the setup method.
  %   C MEX counterpart: mdlCheckParameters
  %
  block.RegBlockMethod('CheckParameters', @CheckPrms);

  %
  % SetInputPortSamplingMode:
  %   Functionality    : Check and set input and output port 
  %                      attributes and specify whether the port is operating 
  %                      in sample-based or frame-based mode
  %   C MEX counterpart: mdlSetInputPortFrameData.
  %   (The DSP System Toolbox is required to set a port as frame-based)
  %
  block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
  
  %
  % SetInputPortDimensions:
  %   Functionality    : Check and set the input and optionally the output
  %                      port dimensions.
  %   C MEX counterpart: mdlSetInputPortDimensionInfo
  %
  block.RegBlockMethod('SetInputPortDimensions', @SetInpPortDims);

  %
  % SetOutputPortDimensions:
  %   Functionality    : Check and set the output and optionally the input
  %                      port dimensions.
  %   C MEX counterpart: mdlSetOutputPortDimensionInfo
  %
  block.RegBlockMethod('SetOutputPortDimensions', @SetOutPortDims);
  
  %
  % SetInputPortDatatype:
  %   Functionality    : Check and set the input and optionally the output
  %                      port datatypes.
  %   C MEX counterpart: mdlSetInputPortDataType
  %
  block.RegBlockMethod('SetInputPortDataType', @SetInpPortDataType);
  
  %
  % SetOutputPortDatatype:
  %   Functionality    : Check and set the output and optionally the input
  %                      port datatypes.
  %   C MEX counterpart: mdlSetOutputPortDataType
  %
  block.RegBlockMethod('SetOutputPortDataType', @SetOutPortDataType);
  
  %
  % SetInputPortComplexSignal:
  %   Functionality    : Check and set the input and optionally the output
  %                      port complexity attributes.
  %   C MEX counterpart: mdlSetInputPortComplexSignal
  %
  block.RegBlockMethod('SetInputPortComplexSignal', @SetInpPortComplexSig);
  
  %
  % SetOutputPortComplexSignal:
  %   Functionality    : Check and set the output and optionally the input
  %                      port complexity attributes.
  %   C MEX counterpart: mdlSetOutputPortComplexSignal
  %
  block.RegBlockMethod('SetOutputPortComplexSignal', @SetOutPortComplexSig);
  
  %
  % PostPropagationSetup:
  %   Functionality    : Set up the work areas and the state variables. You can
  %                      also register run-time methods here.
  %   C MEX counterpart: mdlSetWorkWidths
  %
  block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);

  % -----------------------------------------------------------------
  % Register methods called at run-time
  % -----------------------------------------------------------------
  
  % 
  % ProcessParameters:
  %   Functionality    : Call to allow an update of run-time parameters.
  %   C MEX counterpart: mdlProcessParameters
  %  
  block.RegBlockMethod('ProcessParameters', @ProcessPrms);

  % 
  % InitializeConditions:
  %   Functionality    : Call to initialize the state and the work
  %                      area values.
  %   C MEX counterpart: mdlInitializeConditions
  % 
  block.RegBlockMethod('InitializeConditions', @InitializeConditions);
  
  % 
  % Start:
  %   Functionality    : Call to initialize the state and the work
  %                      area values.
  %   C MEX counterpart: mdlStart
  %
  block.RegBlockMethod('Start', @Start);

  % 
  % Outputs:
  %   Functionality    : Call to generate the block outputs during a
  %                      simulation step.
  %   C MEX counterpart: mdlOutputs
  %
  block.RegBlockMethod('Outputs', @Outputs);

  % 
  % Update:
  %   Functionality    : Call to update the discrete states
  %                      during a simulation step.
  %   C MEX counterpart: mdlUpdate
  %
  block.RegBlockMethod('Update', @Update);

  % 
  % Derivatives:
  %   Functionality    : Call to update the derivatives of the
  %                      continuous states during a simulation step.
  %   C MEX counterpart: mdlDerivatives
  %
  % not required so remove current implementatin because it is wrong 
  block.RegBlockMethod('Derivatives', @Derivatives);
  
  % 
  % Projection:
  %   Functionality    : Call to update the projections during a
  %                      simulation step.
  %   C MEX counterpart: mdlProjections
  %
  block.RegBlockMethod('Projection', @Projection);
  
  % 
  % SimStatusChange:
  %   Functionality    : Call when simulation enters pause mode
  %                      or leaves pause mode.
  %   C MEX counterpart: mdlSimStatusChange
  %
  block.RegBlockMethod('SimStatusChange', @SimStatusChange);
  
  % 
  % Terminate:
  %   Functionality    : Call at the end of a simulation for cleanup.
  %   C MEX counterpart: mdlTerminate
  %
  block.RegBlockMethod('Terminate', @Terminate);

  %
  % GetOperatingPoint:
  %   Functionality    : Return the operating point of the block.
  %   C MEX counterpart: mdlGetOperatingPoint
  %
  block.RegBlockMethod('GetOperatingPoint', @GetOperatingPoint);
  
  %
  % SetOperatingPoint:
  %   Functionality    : Set the operating point data of the block using
  %                       from the given value.
  %   C MEX counterpart: mdlSetOperatingPoint
  %
  block.RegBlockMethod('SetOperatingPoint', @SetOperatingPoint);

  % -----------------------------------------------------------------
  % Register the methods called during code generation.
  % -----------------------------------------------------------------
  
  %
  % WriteRTW:
  %   Functionality    : Write specific information to model.rtw file.
  %   C MEX counterpart: mdlRTW
  %
  block.RegBlockMethod('WriteRTW', @WriteRTW);
%endfunction

% -------------------------------------------------------------------
% The local functions below are provided to illustrate how you may implement
% the various block methods listed above.
% -------------------------------------------------------------------

function CheckPrms(block)
  
  %a = block.DialogPrm(3).Data;
  %if ~isa(a, 'int32')
  %  me = MSLException(block.BlockHandle, message('Simulink:blocks:invalidParameter'));
  %  throw(me);
  %end
  %a = block.DialogPrm(4).Data;
  %if ~isa(a, 'int32')
  %  me = MSLException(block.BlockHandle, message('Simulink:blocks:invalidParameter'));
  %  throw(me);
  %end
  
%endfunction

function ProcessPrms(block)

  block.AutoUpdateRuntimePrms;
 
%endfunction

function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = fd;
  %block.OutputPort(1).SamplingMode  = fd;
  for outI = 1:1:block.NumOutputPorts
    block.OutputPort(outI).SamplingMode  = fd;
  end
  
%endfunction

function SetInpPortDims(block, idx, di)
  
  block.InputPort(idx).Dimensions = di;
  %block.OutputPort(1).Dimensions  = di;
  for outI = 1:1:block.NumOutputPorts
    block.OutputPort(outI).Dimensions  = di;
  end

%endfunction

function SetOutPortDims(block, idx, di)
  
  block.OutputPort(idx).Dimensions = di;
  %block.InputPort(1).Dimensions    = di;
  for inI = 1:1:block.NumInputPorts
    block.InputPort(inI).Dimensions    = di;
  end

%endfunction

function SetInpPortDataType(block, idx, dt)
  
  block.InputPort(idx).DataTypeID = dt;
  for outI = 1:1:block.NumOutputPorts
    block.OutputPort(outI).DataTypeID = dt;
  end
  %block.OutputPort(1).DataTypeID  = dt;

%endfunction
  
function SetOutPortDataType(block, idx, dt)

  block.OutputPort(idx).DataTypeID  = dt;
  %block.InputPort(1).DataTypeID     = dt;
  for inI = 1:1:block.NumInputPorts
    block.InputPort(inI).DataTypeID     = dt;
  end

%endfunction  

function SetInpPortComplexSig(block, idx, c)
  
  block.InputPort(idx).Complexity = c;
  block.OutputPort(1).Complexity  = c;
  for outI = 1:1:block.NumOutputPorts
    block.OutputPort(outI).Complexity  = c;
  end

%endfunction 
  
function SetOutPortComplexSig(block, idx, c)

  block.OutputPort(idx).Complexity = c;
  block.InputPort(1).Complexity    = c;
  for inI = 1:1:block.NumInputPorts
    block.InputPort(inI).Complexity    = c;
  end

%endfunction 
    
function DoPostPropSetup(block)
  % uncertain if the XyceID and two sim time variables 
  % should be stored in Dwork o
  %display('In DoPostPropSetup')
  block.NumDworks = 5;

  block.Dwork(1).Name            = 'XyceID';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 6;      % int32
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;
  
  block.Dwork(2).Name            = 'XyceSimTimeRequested';
  block.Dwork(2).Dimensions      = 1;
  block.Dwork(2).DatatypeID      = 0;      % double
  block.Dwork(2).Complexity      = 'Real'; % real
  block.Dwork(2).UsedAsDiscState = true;
  
  block.Dwork(3).Name            = 'XyceSimTimeActual';
  block.Dwork(3).Dimensions      = 1;
  block.Dwork(3).DatatypeID      = 0;      % double
  block.Dwork(3).Complexity      = 'Real'; % real
  block.Dwork(3).UsedAsDiscState = true;
  
  block.Dwork(4).Name            = 'XyceSimTimeFinal';
  block.Dwork(4).Dimensions      = 1;
  block.Dwork(4).DatatypeID      = 0;      % double
  block.Dwork(4).Complexity      = 'Real'; % real
  block.Dwork(4).UsedAsDiscState = true;
  
  block.Dwork(5).Name            = 'XyceCirHasADC';
  block.Dwork(5).Dimensions      = 1;
  block.Dwork(5).DatatypeID      = 0;      % double
  block.Dwork(5).Complexity      = 'Real'; % real
  block.Dwork(5).UsedAsDiscState = true;
  
  % Register all tunable parameters as runtime parameters.
  block.AutoRegRuntimePrms;

%endfunction

function InitializeConditions(block)

block.ContStates.Data = 1;

%endfunction

function Start(block)
  % this function is called at the start of a simulation. 
  % we should connect with a Xyce instance 
  % store the Xyce session ID 
  % Have Xyce process in the input netlist and get ready to start.
  display('in Start...')
  if(0)
    % REST access method
    jsarg=jsonencode(' ');
    status = webwritenoproxy("http://localhost:5000/xyce_open", jsarg);
    xyceID = status.Body.Data.uuid;
    %display(['In Start and xyce id is ' status.Body.Data.uuid ] );
    block.Dwork(1).Data = uint32(str2num(xyceID));
    s.uuid=xyceID;
  else
    % python access method
    %pyXyceObj = py.xyce_interface.xyce_interface();
    %display(pyXyceObj);
    %block.Dwork(1).Data = pyXyceObj.xycePtr;
    %block.Dwork(5).Data = pyXyceObj.lib;
    % check if there is a xyce container in the python environment
    if( ~pyrun("a = 'xyceObjects' in locals()", "a"))
      pyrun("xyceObjects = []");
    end
    xyceIdx = pyrun("from xyce_interface import xyce_interface; xyceObjects.append(xyce_interface()); idx = len(xyceObjects)", "idx");
    block.Dwork(1).Data = xyceIdx.int32 - 1;
  end
  
  if(0)
    % REST access method 
    s.simfile=block.DialogPrm(1).Data;
    s.simdir=block.DialogPrm(2).Data;
    jsarg=jsonencode(s);
    status = webwritenoproxy("http://localhost:5000/xyce_initialize", jsarg);
    % should check return code... testCase.verifyEqual(status.StatusCode, matlab.net.http.StatusCode(200));
    % sim time last requested by simulink 
  else
    % python access method
    % use file name for input to Xyces initialize function.
    argv=py.list({block.DialogPrm(1).Data});
    %pyXyceObj.xycePtr = block.Dwork(1).Data;
    %pyXyceObj.lib = block.Dwork(5).Data;
    xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
    status = xyceObj.initialize(argv);
  end 
  
  block.Dwork(2).Data = 0;
  % store -1 in XyceSimTimeActual to indicate that the simulator has not yet 
  % started a transient run nor done the DC op. 
  block.Dwork(3).Data = -1;
  
  if(0)
    % REST access method 
    status = webwritenoproxy("http://localhost:5000/xyce_getfinaltime", jsarg);
    finalTime = status.Body.Data.time;
    block.Dwork(4).Data = finalTime;
  else
    % python access method
    %pyXyceObj = block.Dwork(1).Datal;
    xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
    finalTime = xyceObj.getFinalTime();
    block.Dwork(4).Data = finalTime;
  end 
  
  %display(block.Dwork(4).Data);
  
  % loop over output ports and check for ADC devices as these take a 
  % differnt type of call to gather output from Xyce
  % like the input ports, block.DialogPrm(6).Data is a cell array of the output ports
  % use Xyce results to populate Simulink output.
  
  % assume no ADCs listed in the output
  block.Dwork(5).Data = 0;
  [tr,tc] = size(block.DialogPrm(6).Data );
  for( row = 1:1:tr)
    if( (block.DialogPrm(6).Data{row,1} ~= '-') && (str2num(block.DialogPrm(6).Data{row,1}) <= block.NumInputPorts ) )
      deviceName = block.DialogPrm(6).Data{row,2};
      if ( contains( deviceName, 'ADC'))
        block.Dwork(5).Data = 1;
      end
    end
  end
%endfunction

function WriteRTW(block)
  
   block.WriteRTWParam('matrix', 'M',    [1 2; 3 4]);
   block.WriteRTWParam('string', 'Mode', 'Auto');
   
%endfunction

function Outputs(block)
  % this is called when simulink wants an update of outputs
  % input type is Simulink.MSFcnRunTimeBlock 
  %display(block.CurrentTime);
  
  % copy Simulink inputs to Xyce 
  % block.DialogPrm(5).Data is a cell array of the input data
  % for example:
  %  {'-'}    {'TEMP'            }
  %  {'1'}    {'YDAC!DAC_DRIVER1'}
  %  {'2'}    {'YDAC!DAC_DRIVER2'}
  %  {'-'}    {'YDAC!DAC_DRIVER3'}
  % 
  % the first colum is the simulink port and the second is the Xyce name 
  % ports with "-" are not used 
  
  if (block.CurrentTime < block.Dwork(4).Data)
    [tr,tc] = size(block.DialogPrm(5).Data );
    for( row = 1:1:tr)
      if( (block.DialogPrm(5).Data{row,1} ~= '-') && (str2num(block.DialogPrm(5).Data{row,1}) <= block.NumInputPorts ) )
        %display( block.DialogPrm(5).Data{row,1})
        portID = str2num(block.DialogPrm(5).Data{row,1});
        %display( block.DialogPrm(5).Data{row,2})
        deviceName = block.DialogPrm(5).Data{row,2};
        
        % need to make this more general because setting a DAC value and a general
        % circuit parameter value use different methods.  
        if ( contains( deviceName, 'DAC'))
          % need at least two points here as it should be an array. But do the time points need to be unique?
          timeArray = [block.CurrentTime ; block.CurrentTime ];
          dacVData = [block.InputPort(portID).Data; block.InputPort(portID).Data];
          if(0)
            % REST access method
            % get the xyceID
            s.uuid=num2str(block.Dwork(1).Data);
            s.devname = deviceName;
            s.timearray = timeArray;
            s.voltarray = dacVData;
            jsarg=jsonencode(s);
            status = webwritenoproxy("http://localhost:5000/xyce_updatetimevoltagepairs", jsarg);
          else
            % python access method 
            xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
            status = xyceObj.updateTimeVoltagePairs( deviceName, timeArray, dacVData);
          end
        else
          % not a DAC so use the more generic setcircuit parameter method
          if(0)
            % REST access method
            s.uuid=num2str(block.Dwork(1).Data);
            s.paramname=deviceName;
            s.paramval=block.InputPort(portID).Data;
            jsarg=jsonencode(s);
            status = webwritenoproxy("http://localhost:5000/xyce_setcircuitparameter", jsarg);
          else
            % python access method 
            xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
            status = xyceObj.setCircuitParameter(deviceName, block.InputPort(portID).Data);
          end
          
        end 
      end
    end
  end
      
      
  % check that CurrentTime is past where the simulator is but not further than 
  % where the simulator is configured to go.
  %
  % note Dwork(3) = XyceSimTimeActual
  %      Dwork(3) = XyceSimTimeFinal
  %
  % should probably warn user if or when simulink time exceedes XyceSimTime
  % or fix Xyce so one can set final sim time from the API
  while ((block.Dwork(3).Data < block.CurrentTime) && (block.CurrentTime < block.Dwork(4).Data) )
    % try to advance time to block.CurrentTime 
    % get the Xyce uuid
    if(0)
      % REST access method 
      s.uuid = num2str(block.Dwork(1).Data);
      s.simtime = block.CurrentTime;
      jsarg=jsonencode(s);
      status = webwritenoproxy("http://localhost:5000/xyce_simulateuntil", jsarg);
      %
      % still need to handle case where Xyce fails on a timestep and cannot continue
      % 
      %display(status);
      block.Dwork(3).Data = status.Body.Data.simulatedTime;
      %testCase.verifyEqual(retcode.simulatedTime, (i * deltaTime), "AbsTol", 1e-6);
    else 
      % python access method 
      xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
      % this returns a python tuple
      ptAns = xyceObj.simulateUntil( block.CurrentTime );
      % need to convert the tuple to a cell array in Matlab
      cptAns = cell(ptAns);
      % now one can access the data returned from python
      block.Dwork(3).Data = cptAns{2};
    end 
  end
  
  % like the input ports, block.DialogPrm(6).Data is a cell array of the output ports
  % use Xyce results to populate Simulink output.
  if (block.CurrentTime < block.Dwork(4).Data)
    % all of the ADC values in a circuit are obtained with one call.
    % if the user is asking for an ADC output then grab all of the ADC data
    adcCircuitData = struct;
    if(block.Dwork(5).Data == 1)
      if(0)
        % REST access method
        s.uuid = num2str(block.Dwork(1).Data);
        jsarg=jsonencode(s);
        status = webwritenoproxy("http://localhost:5000/xyce_gettimevoltagepairs", jsarg);
        % status.Body.Data is a struct with the following info
        % ADCnames {} a cell array
        % numADCnames : number of returned ADC names
        % numPointsInArray : an array with the number of returned data points numADCs by num time points
        % timeArray : an array of time values returned 
        % voltageArray  : voltage levels of the ADCs 
        adcCircuitData = status.Body.Data;
      else 
        % python access method 
        xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
        % (status, ADCnames, numADCnames, numPointsArray, timeArray, voltageArray)  = xyceObj.getTimeVoltagePairsADCLimitData()
        ptAns = xyceObj.getTimeVoltagePairsADCLimitData();
        cptAns = cell(ptAns);
        adcCircuitData.ADCnames = string(cptAns{2});
        adcCircuitData.numADCnames = int32(cptAns{3});
        adcCircuitData.numPointsInArray = int32(cptAns{4});
        adcCircuitData.timeArray = cptAns{5};
        adcCircuitData.voltageArray = cptAns{6};
      end 
    end
    
    [tr,tc] = size(block.DialogPrm(6).Data );
    for( row = 1:1:tr)
      if( (block.DialogPrm(6).Data{row,1} ~= '-') && (str2num(block.DialogPrm(6).Data{row,1}) <= block.NumOutputPorts ) )
        portNum = str2num(block.DialogPrm(6).Data{row,1});
        deviceName = block.DialogPrm(6).Data{row,2};
        if ( contains( deviceName, 'ADC'))
            for j = 1:adcCircuitData.numADCnames
              if (adcCircuitData.ADCnames{j} == deviceName)
                numPoints = int32(adcCircuitData.numPointsInArray(j)); 
                block.OutputPort(portNum).Data = adcCircuitData.voltageArray{j}{numPoints};
              end 
            end
        else
          if(0)
            % REST access method 
            % get the Xyce uuid
            s.uuid = num2str(block.Dwork(1).Data);
            s.paramname =deviceName;
            jsarg=jsonencode(s);
            status = webwritenoproxy("http://localhost:5000/xyce_getcircuitvalue", jsarg);
            %display( status.Body. );
            block.OutputPort(portNum).Data = status.Body.Data.value;
          else 
            % python access method 
            xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
            block.OutputPort(portNum).Data = xyceObj.getCircuitValue(deviceName);
          end
        end
      end
    end
  end  
  %block.OutputPort(1).Data = 1.5*block.InputPort(1).Data;
  
%endfunction

function Update(block)
  % less clear about when this is called.  According to the documentation
  % it is called at major time steps.
  %display('In Update...')
  block.Dwork(2).Data = block.InputPort(1).Data;
  
%endfunction

function Derivatives(block)

  block.Derivatives.Data = block.ContStates.Data;

%endfunction

function Projection(block)

  states = block.ContStates.Data;
  block.ContStates.Data = states+eps; 

%endfunction

function SimStatusChange(block, s)
  
  block.Dwork(2).Data = block.Dwork(2).Data+1;    

  if s == 0
    disp('Pause in simulation.');
  elseif s == 1
    disp('Resume simulation.');
  end
  
%endfunction
    
function Terminate(block)
  %disp(['Terminating the block with handle ' num2str(block.BlockHandle) '.']);
  %display(['In terminate. ' num2str(block.Dwork(1).Data) ]);
  if(0)
    % REST access method
    s.uuid = num2str(block.Dwork(1).Data);
    %display(['In terminate pt2. ' s.uuid ]);
    jsarg=jsonencode(s);
      
    % close this xyce object
    status = webwritenoproxy("http://localhost:5000/xyce_close", jsarg);
    % display error message if this fails. display(status)
  else 
    % python access method
    xyceObj = pyrun("x = xyceObjects[ int(idx)]", "x", idx = block.Dwork(1).Data);
    %close the object
    xyceObj.close();
    %remove it from the array
    pyrun("del xyceObjects[int(idx)]", idx=block.Dwork(1).Data);
  end
%endfunction
 
function operPointData = GetOperatingPoint(block)
  % package the Dwork data as the entire operating point of this block
  operPointData = block.Dwork(1).Data;

%endfunction

function SetOperatingPoint(block, operPointData)
  % the operating point of this block is the Dwork data (this method 
  % typically performs the inverse actions of the method GetOperatingPoint)
  block.Dwork(1).Data = operPointData;

%endfunction
