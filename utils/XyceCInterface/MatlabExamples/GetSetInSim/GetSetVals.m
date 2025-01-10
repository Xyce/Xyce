



py.importlib.import_module('xyce_interface');

pyXyceObj = py.xyce_interface.xyce_interface();


argv=py.list({"TestNetlist2.cir"});

pyXyceObj.initialize(argv);

r1res = pyXyceObj.getCircuitValue('R1:R');
display("==> R1 R=" + r1res);

cirTemp = pyXyceObj.getCircuitValue("TEMP");
display("==> Circuit Temp = " + cirTemp);


numSteps=100;
deltaTime = 1.0 / numSteps;
for i = 0:numSteps
  if( i==(numSteps/2)) 
    % in the middle of the run change the resistance on resistor R1
    setParamResult = pyXyceObj.setCircuitParameter( "R1:R", 1000.0);
    setParamResult = pyXyceObj.setCircuitParameter( "TEMP", -50.0);
  end

  if( i==((numSteps/2)+1))
    % in the middle of the run change the resistance on resistor R1
    circuitParamValue = 0.0;
    circuitParamValue = pyXyceObj.getCircuitValue( "R1:R");
    display( "==> R1 R = " + circuitParamValue);
    circuitParamValue = pyXyceObj.getCircuitValue( "TEMP");
    display( "==> Circuit temperature is now = " + circuitParamValue);
  end
 
  display("==> Taking step " + i + " of " + numSteps);
  requested_time = 0.0 + (i+1) * deltaTime;
  actual_time = 0.0;
  result  = pyXyceObj.simulateUntil( requested_time );
end
% get some measure results at the end of the run

maxV1=pyXyceObj.getCircuitValue( "MAXV1");
minV1=pyXyceObj.getCircuitValue( "MINV1");
maxV2=pyXyceObj.getCircuitValue( "MAXV2");
minV2=pyXyceObj.getCircuitValue( "MINV2");

display("==> Measure results from Xyce V1 = " + maxV1 + ", " + minV1 + " V2 = " + maxV2 + "," + minV2);


pyXyceObj.runSimulation();

pyXyceObj.close();