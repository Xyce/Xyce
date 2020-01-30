//-------------------------------------------------------------------------
//   Copyright 2002-2020 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
//
//
//
//-------------------------------------------------------------------------


#ifndef _CHARON_CIRCUIT_INTERFACE_H_
#define _CHARON_CIRCUIT_INTERFACE_H_

#include <map>
#include <vector>
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParamList;
}

namespace charon {

  namespace sc {

    /**
     * @brief - Provides the interface for external circuit codes to call charon to take a transient step.
     *
     * This is a singleton so we can avoid punching layers into xyce
     * to pass this in.
     */
    class CircuitInterface {
      
    public:
      
      //! Destructor.
      ~CircuitInterface();

      //! Returns an instance of this object. 
      static CircuitInterface& getInstance();

      //! Calls the CharonClient::Run_Step_CircuitSimulation routine to start a charon solve.
      bool takeStep(
            const Teuchos::RCP<Teuchos::ParameterList>& inputList,
	    const std::map<std::string, double>& inputMap,
	    const Teuchos::RCP<Teuchos::ParameterList>& outputList,
	    std::vector<double>& outputVector,
	    std::vector< std::vector<double> >& outputJacobian);

      //! Tells Charon to accept the time step and write output.
      void acceptTimeStep(bool& is_active);

    private:
      
      /*! Tells Charon to reset the time integrator and that we will be making another solve at the current time step.  This should be called after every solve of charon unless we accept the step, otherwise, the transient history will be corrupted.  When Charon attempts a step, it pushes the new solution and time step size onto the stored history stack.  If the step fails we need to remove it from the stack.
       */
      void declineTimeStep();

      //! Singleton - disallow constructor.
      CircuitInterface();

      //! Singleton - disallow copy constructor.
      CircuitInterface(const CircuitInterface& source) {};

    protected:

      //! Used to determine if the next step should be reset. 
      bool accepted_step_;

    };
    
  } // END "namespace sc"
  
} // END "namespace charon"

#endif 
