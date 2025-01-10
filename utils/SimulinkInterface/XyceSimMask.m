classdef XyceSimMask

  methods(Static)

    
    function MaskInitialization(maskInitContext)
      % Following properties of 'maskInitContext' are available to use:
      %  - BlockHandle 
      %  - MaskObject 
      %  - MaskWorkspace: Use get/set APIs to work with mask workspace.
      blkHandle = maskInitContext.BlockHandle;       % Block Handle of this block
      maskObj = maskInitContext.MaskObject;          % Mask object of this masked block
      maskWorkspace = maskInitContext.MaskWorkspace; % Use this to work with mask workspace
      %disp('Mask initialization triggered');
      % get the number of inputs and outputs as saved in the I/O selectors
      numInputsSpinnerValue = str2num(get_param(gcb, 'NumberOfInputs'));
      maskWorkspace.set('NumInputPorts', numInputsSpinnerValue);
      numOutputsSpinnerValue = str2num(get_param(gcb, 'NumberOfOutputs'));
      maskWorkspace.set('NumberOfOutputs', numOutputsSpinnerValue);
      
      %end            
    end
    
    % Use the code browser on the left to add the callbacks.
    function SelectXyceInputFileButton(callbackContext)
      [aFile, aPath] = uigetfile('*');
      %display(aFile);
      %display(aPath);
      if( aFile ~= 0 )
        set_param( gcb, 'XyceInputFileName', aFile );
        % picked a new file.  It's a good idea to clear the
        % inputs and outputs tables
        maskObj = get_param( gcb, 'MaskObject');
        tableControl = maskObj.getDialogControl('InputNames');
        for j = tableControl.getNumberOfRows():-1:1
          tableControl.removeRow(j);
        end
        %tableControl.addRow('-','TEMP');

      end
      if( aPath ~= 0 )
        set_param( gcb, 'WorkingDirectory', aPath);
      end
    end

    function ResetForm(callbackContext)
      % reset all inputs on form to defaults
      maskObj = get_param( gcb, 'MaskObject');
      set_param( gcb, 'XyceInputFileName', '<name of circuit file>' );
      set_param( gcb, 'WorkingDirectory', '<path to circuit file>');
      
      tableControl = maskObj.getDialogControl('InputNames');
      for j = tableControl.getNumberOfRows():-1:1
        tableControl.removeRow(j);
      end
      %tableControl.addRow('-','TEMP');
      %set_param(gcb, 'NumberOfInputs', num2str(1));

      tableControl = maskObj.getDialogControl('OutputNames');
      for j = tableControl.getNumberOfRows():-1:1
        tableControl.removeRow(j);
      end
      %tableControl.addRow('-','TEMP');
      %set_param(gcb, 'NumberOfOutputs', num2str(1));
    end
   
    function ScanXyceInputFileButton(callbackContext)
      maskObj = get_param( gcb, 'MaskObject');
      inputFileName  = get_param( gcb, 'XyceInputFileName');
      if( 0 )
        %potentialInputNames = getDACNamesFromCircuit(inputFileName);
        s.stats=0;
        jsarg=jsonencode(s);
        status = webwritenoproxy("http://localhost:5000/xyce_open", jsarg);
        xyceID = status.Body.Data.uuid;
      
        % call initialize 
        s2.uuid=xyceID;
        s2.simfile=inputFileName;
        jsarg2=jsonencode(s2);
        status = webwritenoproxy("http://localhost:5000/xyce_initialize", jsarg2);   
        % call xyce_getdacnames 
        % now try running a circuit.
        s3.uuid=xyceID;
        jsarg3=jsonencode(s3);
        status = webwritenoproxy("http://localhost:5000/xyce_getdacnames", jsarg3);
        dacNames = [];
        if( status.Body.Data.status == 1.0)
          dacNames = status.Body.Data.dacNames;
        end
        status = webwritenoproxy("http://localhost:5000/xyce_getadcmap", jsarg3);
        adcNames = [];
        if( status.Body.Data.status == 1.0)
          adcNames = status.Body.Data.ADCnames;
        end
      
        % close this xyce object
        status = webwritenoproxy("http://localhost:5000/xyce_close", jsarg3);
      end 
      
      % user the python interface connection 
      pyXyceObj = py.xyce_interface.xyce_interface();
      % should check that underlying Xyce pointer is not None at this time.
      argv=py.list({inputFileName});
      status = pyXyceObj.initialize(argv);
      pStatus = pyXyceObj.getDACDeviceNames();
      dacNames = cell(pStatus{2});
      pStatus = pyXyceObj.getADCMap();
      adcNames = cell(pStatus{2});
      pStatus = pyXyceObj.close();
      
      % this is a controller that lets us add to the table
      tableControl = maskObj.getDialogControl('InputNames');
      numDacs = size(dacNames,2);
      if( numDacs > 0)
        for r = numDacs:-1:1
          aDacName = string(dacNames{r});
          alreadyFound = 0;
          for j = 1:tableControl.getNumberOfRows()
            cellValue = tableControl.getValue([j 2]);
            if strcmp(aDacName, cellValue )
              alreadyFound = 1;
              break;
            end
          end
          if( ~alreadyFound)
            tableControl.addRow('-', aDacName);
          end
        end
      else 
        for j = tableControl.getNumberOfRows():-1:1
            cellValue = tableControl.getValue([j 2]);
            if ~strcmp('TEMP', cellValue )
              tableControl.removeRow(j);
            end
        end
      end
     
      % this is a controller that lets us add to the table
      tableControl = maskObj.getDialogControl('OutputNames');
      numAdcs = size(adcNames,2)
      if( numAdcs > 0)
        for r = numAdcs:-1:1
          aName = string(adcNames{r});
          alreadyFound = 0;
          for j = 1:tableControl.getNumberOfRows()
            cellValue = tableControl.getValue([j 2]);
            if strcmp(aName, cellValue )
              alreadyFound = 1;
              break;
            end
          end
          if( ~alreadyFound)
            tableControl.addRow('-', aName);
          end
        end
      else 
        for j = tableControl.getNumberOfRows():-1:1
            cellValue = tableControl.getValue([j 2]);
            if ~strcmp('TEMP', cellValue )
              tableControl.removeRow(j);
            end
        end
      end

    end

   
   
    function NumberOfInputs(callbackContext)
      % when the number of inputs is changed, update
      % the potential port nubers in the InputNames Table
      maskObj = get_param( gcb, 'MaskObject');
      numInputsSpinnerValue = str2num(get_param(gcb, 'NumberOfInputs'));
      %set_param( gcb, 'NumberOfInputPorts', get_param(gcb, 'NumberOfInputs' );
      %set_param(gcb, 'XyceNumberOfInputPorts', numInputsSpinnerValue);
      
      possibleInputs = {'-' '1'};
      if (numInputsSpinnerValue > 1)
        for i = 1:1:(numInputsSpinnerValue)
            possibleInputs{i+1} = int2str(i);
        end
      end
      tableControl = maskObj.getDialogControl('InputNames');
      if( tableControl.getNumberOfColumns > 1)
        portNumberColumn = tableControl.getColumn(1);
        portNumberColumn.TypeOptions = possibleInputs;
        if (tableControl.getNumberOfRows() < numInputsSpinnerValue)
          for j = tableControl.getNumberOfRows():1:(numInputsSpinnerValue)
            tableControl.addRow('-','-');
          end
        elseif (tableControl.getNumberOfRows() > numInputsSpinnerValue)
          for j = tableControl.getNumberOfRows():-1:(numInputsSpinnerValue+1)
            tableControl.removeRow(j);
          end
        end
      end 
    end

    function NumberOfOutputs(callbackContext)
      % when the number of outputs is changed, update
      % the potential port nubers in the InputNames Table
      maskObj = get_param( gcb, 'MaskObject');
      numOutputsSpinnerValue = str2num(get_param(gcb, 'NumberOfOutputs'));
      
      %set_param(gcb, 'XyceNumberOfOutputPorts', numOutputsSpinnerValue);
      possibleOutputs = {'-' '1'};
      if (numOutputsSpinnerValue > 1)
        for i = 1:1:(numOutputsSpinnerValue)
            possibleOutputs{i+1} = int2str(i);
        end
      end
      tableControl = maskObj.getDialogControl('OutputNames');
      if( tableControl.getNumberOfColumns > 1)
        portNumberColumn = tableControl.getColumn(1);
        portNumberColumn.TypeOptions = possibleOutputs;
        if (tableControl.getNumberOfRows() < numOutputsSpinnerValue)
          for j = tableControl.getNumberOfRows():1:(numOutputsSpinnerValue)
            tableControl.addRow('-','-');
          end
        elseif (tableControl.getNumberOfRows() > numOutputsSpinnerValue)
          for j = tableControl.getNumberOfRows():-1:(numOutputsSpinnerValue+1)
            tableControl.removeRow(j);
          end
        end 
      end
    end

    function InputNames(callbackContext)
      %display('InputNames callback called.');
      maskObj = get_param( gcb, 'MaskObject');
      
      numInputsSpinnerValue = str2num(get_param(gcb, 'NumberOfInputs'));
      tempInputNames = cell(1, numInputsSpinnerValue);
      for i = 1:1:numInputsSpinnerValue
          tempInputNames{i} = '';
      end
      %display(tempInputNames);
      % loop over inputs names setting the port names up for the 
      % first numInputs ports
      tableControl = maskObj.getDialogControl('InputNames');
      numRows = tableControl.getNumberOfRows();
      for i = 1:1:numRows
          portNum = tableControl.getValue([i 1]);
          [val, numCheck] = str2num(portNum);
          if( numCheck == 1)
              % port is a number
              name = tableControl.getValue([i, 2]);
              tempInputNames{val} = name;
          end
          
      end
      %display(tempInputNames);
      xtextArea = get_param( gcb, 'XyceInputPortNames');
      xtextArea.Value = tempInputNames;
      %set_param(gcb, 'XyceInputPortNames', xtextArea);
    end

    function OutputNames(callbackContext)
      %display('OutputNames callback called');
      maskObj = get_param( gcb, 'MaskObject');
      
      numOutputsSpinnerValue = str2num(get_param(gcb, 'NumberOfOutputs'));
      tempOutputNames = cell(1, numOutputsSpinnerValue);
      for i = 1:1:numOutputsSpinnerValue
          tempOutputNames{i} = '';
      end
      %display(tempOutputNames);
      % loop over output names setting the port names up for the 
      % first numOutputs ports
      tableControl = maskObj.getDialogControl('OutputNames');
      numRows = tableControl.getNumberOfRows();
      for i = 1:1:numRows
          portNum = tableControl.getValue([i 1]);
          [val, numCheck] = str2num(portNum);
          if( numCheck == 1)
              % port is a number
              name = tableControl.getValue([i, 2]);
              tempOutputNames{val} = name;
          end
          
      end
      %display(tempOutputNames);
      xtextArea = get_param( gcb, 'XyceOutputPortNames');
      xtextArea.Value = tempOutputNames;
    end

    function LoadCallback(callbackContext)
      disp('LoadCallback')
    end

    function AddInputLineButton(callbackContext)
      maskObj = get_param( gcb, 'MaskObject');
      numInputsSpinnerValue = str2num(get_param(gcb, 'NumberOfInputs'));
      possibleInputs = {'-' '1'};
      if (numInputsSpinnerValue > 1)
        for i = 1:1:(numInputsSpinnerValue)
            possibleInputs{i+1} = int2str(i);
        end
      end
      tableControl = maskObj.getDialogControl('InputNames');
      if( tableControl.getNumberOfColumns > 1)
        %portNumberColumn = tableControl.getColumn(1);
        %portNumberColumn.TypeOptions = possibleOutputs;
        tableControl.addRow('-','-');
      end
    end

    function RemoveInputLineButton(callbackContext)
      maskObj = get_param( gcb, 'MaskObject');
      tableControl = maskObj.getDialogControl('InputNames');
      if( tableControl.getNumberOfRows >= 1)
        numRows = tableControl.getNumberOfRows();
        tableControl.removeRow(numRows);
      end
    end

    function AddOutputLineButton(callbackContext)
      maskObj = get_param( gcb, 'MaskObject');
      numOutputsSpinnerValue = str2num(get_param(gcb, 'NumberOfOutputs'));
      possibleOutputs = {'-' '1'};
      if (numOutputsSpinnerValue > 1)
        for i = 1:1:(numOutputsSpinnerValue)
            possibleOutputs{i+1} = int2str(i);
        end
      end
      tableControl = maskObj.getDialogControl('OutputNames');
      if( tableControl.getNumberOfColumns > 1)
        %portNumberColumn = tableControl.getColumn(1);
        %portNumberColumn.TypeOptions = possibleOutputs;
        tableControl.addRow('-','-');
      end
    end

    function RemoveOutputLineButton(callbackContext)
      maskObj = get_param( gcb, 'MaskObject');
      tableControl = maskObj.getDialogControl('OutputNames');
      if( tableControl.getNumberOfRows >= 1)
        numRows = tableControl.getNumberOfRows();
        tableControl.removeRow(numRows);
      end
      
    end
  end
end