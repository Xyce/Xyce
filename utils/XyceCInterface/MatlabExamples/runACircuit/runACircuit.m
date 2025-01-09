



py.importlib.import_module('xyce_interface');

pyXyceObj = py.xyce_interface.xyce_interface()


argv=py.list({"runACircuit.cir"})

pyXyceObj.initialize(argv)

pyXyceObj.runSimulation()

pyXyceObj.close()