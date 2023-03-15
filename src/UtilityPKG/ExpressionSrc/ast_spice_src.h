//-------------------------------------------------------------------------
//   Copyright 2002-2023 National Technology & Engineering Solutions of
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

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 10/xx/2019
//
//
//
//
//-----------------------------------------------------------------------------

#ifndef _ast_spice_src_h_
#define _ast_spice_src_h_

#define AST_CALL_SUBFUNC(PTR,FUNC,ARG)  if(this->PTR) { this->PTR->FUNC(ARG);  }

#define AST_CALL_SUBOUTPUT(PTR)  if(  !(Teuchos::is_null(PTR))  ) { os << std::setw(indent) << " "; os << #PTR << ": " << std::endl; PTR->output(os,indent+1); }

//-------------------------------------------------------------------------------
// spice pulse  operator
//
// pulse( v1,v2,td,tr,tf,pw,per)
//
template <typename ScalarT>
class spicePulseOp : public astNode<ScalarT>
{
  public:
    spicePulseOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(args), 
      time_(time), 
      v1Given_(false), v2Given_(false), tdGiven_(false),
      trGiven_(false), tfGiven_(false), pwGiven_(false), perGiven_(false),
      bpTol_(0.0), startingTimeStep_(0.0), finalTime_(0.0)
  {
    if (args.size() < 1)
    {
      std::vector<std::string> errStr(1,std::string("AST node (spice_pulse) needs at least 1 argument.  V1 is required for the PULSE source function.")); yyerror(errStr);
    }

    if (args.size() > 7)
    {
      std::vector<std::string> errStr(1,std::string("AST node (spice_pulse) has too many arguments")); yyerror(errStr);
    }


    // ERK. Fix this to set the proper defaults, consistent with N_DEV_PulseData
    //
    // double tstep = solState_.startingTimeStep_;
    // double tstop = solState_.finalTime_;
    // if (!TDgiven)  TD  = 0.0;
    // if (!TRgiven)  TR  = tstep;
    // if (!TFgiven)  TF  = tstep;
    // if (!PWgiven)  PW  = tstop;
    // if (!PERgiven) PER = tstop;
    //
    // At the time of construction, I don't think I have these values yet, so set in the val function.
    //
    std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
    if (child.size() < 7) { child.resize(7); }

    if (args.size() >= 1) { v1Given_=true; } else { child[0] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
    if (args.size() >= 2) { v2Given_=true; } else { child[1] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
    if (args.size() >= 3) { tdGiven_=true; } else { child[2] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
    if (args.size() >= 4) { trGiven_=true; } else { child[3] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
    if (args.size() >= 5) { tfGiven_=true; } else { child[4] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
    if (args.size() >= 6) { pwGiven_=true; } else { child[5] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
    if (args.size() >= 7) { perGiven_=true; } else { child[6] = ( Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)) ); }
  };

    virtual ScalarT val()
    {
      if (!trGiven_) { Teuchos::RCP<numval<ScalarT> > trTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (this->childrenAstNodes_[3]); trTmpOp->number = startingTimeStep_; }
      if (!tfGiven_) { Teuchos::RCP<numval<ScalarT> > tfTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (this->childrenAstNodes_[4]); tfTmpOp->number = startingTimeStep_; }
      if (!pwGiven_) { Teuchos::RCP<numval<ScalarT> > pwTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (this->childrenAstNodes_[5]); pwTmpOp->number = finalTime_; }
      if (!perGiven_) { Teuchos::RCP<numval<ScalarT> > perTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (this->childrenAstNodes_[6]); perTmpOp->number = finalTime_; }

      ScalarT time = std::real(this->time_->val());
      ScalarT V1 = this->childrenAstNodes_[0]->val();
      ScalarT V2 = this->childrenAstNodes_[1]->val();
      ScalarT TD = std::real(this->childrenAstNodes_[2]->val());
      ScalarT TR = std::real(this->childrenAstNodes_[3]->val());
      ScalarT TF = std::real(this->childrenAstNodes_[4]->val());
      ScalarT PW = std::real(this->childrenAstNodes_[5]->val());
      ScalarT PER = std::real(this->childrenAstNodes_[6]->val());

      time -= TD;

      if (std::real(time) > std::real(PER) && std::real(PER) != 0.0)
      {
        // repeating signal - figure out where we are in period
        ScalarT basetime = std::real(PER) * std::floor(std::real(time)/std::real(PER));
        time -= basetime;
      }

      // This section got ugly because of a nasty roundoff bug.
      // Instead of doing "time > X" you need also check that time
      // is not within bptol of X.
      // So the following translation is used:
      // Instead of:                           we do:
      //  time > X                            time>X && fabs(time-x)>bptol
      //  time <= X                           time<X || fabs(time-x)<bptol

      ScalarT SourceValue = 0.0;

      if (std::real(time) <= 0 || (std::real(time) > (std::real(TR) + std::real(PW) + std::real(TF)) &&
            (fabs (std::real(time) - std::real(TR+PW+TF)) > std::real(bpTol_)) ) )
      {
        SourceValue = V1;
      }
      else if ((std::real(time) > std::real(TR) && fabs(std::real(time)-std::real(TR)) > std::real(bpTol_))
          && (std::real(time) < (std::real(TR) + std::real(PW)) || fabs (std::real(time)-std::real(TR+PW))<std::real(bpTol_)) )
      {
        SourceValue = V2;
      }
      else if (std::real(time) > 0 && (std::real(time) < std::real(TR) || fabs(std::real(time)-std::real(TR)) < std::real(bpTol_)))
      {
        if (std::real(TR) != 0.0)
        {
          SourceValue = V1 + (V2 - V1) * (std::real(time)) / std::real(TR);
        }
        else
        {
          SourceValue = V1;
        }
      }
      else
      { // time > (TR + PW) && <= (TR + PW + TF)
        if (std::real(TF) != 0.0)
        {
          SourceValue = V2 + (V1 - V2) * (std::real(time) - (std::real(TR) + std::real(PW))) / std::real(TF);
        }
        else
        {
          SourceValue = V2; //      SourceValue = 0.5 * (V1 + V2);
        }
      }

      return SourceValue;
    }

    virtual ScalarT dx (int i)
    {
      return 0.0;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) 
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, only used for transient. update if this changes.
    virtual bool getIsComplex () { return false; }

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      ScalarT time = std::real(this->time_->val());
      ScalarT V1 = this->childrenAstNodes_[0]->val();
      ScalarT V2 = this->childrenAstNodes_[1]->val();
      ScalarT TD = std::real(this->childrenAstNodes_[2]->val());
      ScalarT TR = std::real(this->childrenAstNodes_[3]->val());
      ScalarT TF = std::real(this->childrenAstNodes_[4]->val());
      ScalarT PW = std::real(this->childrenAstNodes_[5]->val());
      ScalarT PER = std::real(this->childrenAstNodes_[6]->val());

      int currPeriodIndex = 0;
      double basetime = 0.0;
      time -= TD;

      // repeating signal - figure out where we are in period
      if(std::real(time) >= std::real(PER))
      {
        if (std::real(PER) != 0.0)
        {
          currPeriodIndex = (static_cast<int> (floor(std::real(time)/std::real(PER))));
          basetime = std::real(PER) * (static_cast<double> (currPeriodIndex));
          time -= basetime;
        }
      }

      // now that we know which period this is, push_back all breakpoints
      // in this period and the next.  If we are still in the delay, then
      // just use first two periods.

      // current period:
      breakPointTimes.push_back(std::real(basetime+TD));
      breakPointTimes.push_back(std::real(basetime+TD+TR));
      breakPointTimes.push_back(std::real(basetime+TD+TR+PW));
      breakPointTimes.push_back(std::real(basetime+TD+TR+PW+TF));

      if (std::real(PER) != 0.0)
      {
        breakPointTimes.push_back(std::real(basetime+TD+PER));

        // next period:
        breakPointTimes.push_back(std::real(basetime+TD+PER+TR));
        breakPointTimes.push_back(std::real(basetime+TD+PER+TR+PW));
        breakPointTimes.push_back(std::real(basetime+TD+PER+TR+PW+TF));
        breakPointTimes.push_back(std::real(basetime+TD+PER+PER));
      }

      return true;
    }

    virtual void setBreakPointTol(double tol) { bpTol_ = tol; }
    virtual void setStartingTimeStep(double timeStep) { startingTimeStep_ = timeStep; }
    virtual void setFinalTime(double finalTime) { finalTime_ = finalTime; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice pulse operator id = " << this->id_ << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[0])
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[1])
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[2])
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[3])
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[4])
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[5])
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[6])
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "spice pulse operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_pulse codeGen function is not implemented yet" <<std::endl;
    }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool srcType() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      this->thisAstNode_ = thisAst_;
      Teuchos::RCP<spicePulseOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<spicePulseOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
                                   //
      this->childrenAstNodes_[0]->accept(visitor,this->childrenAstNodes_[0]);
      this->childrenAstNodes_[1]->accept(visitor,this->childrenAstNodes_[1]);
      this->childrenAstNodes_[2]->accept(visitor,this->childrenAstNodes_[2]);
      this->childrenAstNodes_[3]->accept(visitor,this->childrenAstNodes_[3]);
      this->childrenAstNodes_[4]->accept(visitor,this->childrenAstNodes_[4]);
      this->childrenAstNodes_[5]->accept(visitor,this->childrenAstNodes_[5]);
      this->childrenAstNodes_[6]->accept(visitor,this->childrenAstNodes_[6]);

      time_->accept(visitor,time_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > time_;
    bool v1Given_, v2Given_, tdGiven_, trGiven_, tfGiven_, pwGiven_, perGiven_;
    double bpTol_;
    double startingTimeStep_;
    double finalTime_;
};

//-------------------------------------------------------------------------------
// spice sin operator
//
// sin (v0,va,freq,td,theta,phase)
//
template <typename ScalarT>
class spiceSinOp : public astNode<ScalarT>
{
  public:
    spiceSinOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(args), 
      time_(time),
      v0Given_(false), vaGiven_(false), freqGiven_(false), tdGiven_(false), thetaGiven_(false), phaseGiven_(false),
      finalTime_(0.0)
    {
      if (args.size() < 3)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sin) needs at least 3 argument.  V0, VA and FREQ are required for the SIN source function.")); yyerror(errStr);
      }

      if (args.size() > 6)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sin) has too many arguments")); yyerror(errStr);
      }
   
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      if (child.size() < 6) { child.resize(6); }

      if (args.size() >= 1) { v0Given_=true;    } else { child[0] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 2) { vaGiven_=true;    } else { child[1] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 3) { freqGiven_=true;  } else { child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 4) { tdGiven_=true;    } else { child[3] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 5) { thetaGiven_=true; } else { child[4] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 6) { phaseGiven_=true; } else { child[5] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    };

    virtual ScalarT val()
    {
      // defaults are as follows: (have to be set here b/c in constructor not sure if we know the final time yet
      //
      //   double tstop = solState_.finalTime_;
      //  if (!FREQgiven)  FREQ  = 1.0/tstop;

      if (!freqGiven_ && finalTime_ != 0.0)  
      {
        Teuchos::RCP<numval<ScalarT> > freqTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (this->childrenAstNodes_[2]);
        freqTmpOp->number = 1.0/finalTime_;
      }

      ScalarT time = std::real(this->time_->val());
      ScalarT tdVal = this->childrenAstNodes_[3]->val();
      time -= std::real(tdVal);
      double mpi = M_PI;
      ScalarT SourceValue = 0.0;

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      ScalarT v0Val    = child[0]->val();
      ScalarT vaVal    = child[1]->val();
      ScalarT freqVal  = child[2]->val();
      //ScalarT tdVal    = child[3]->val();
      ScalarT thetaVal = child[4]->val();
      ScalarT phaseVal = child[5]->val();

      if (std::real(time) <= 0)
      {
        SourceValue = (v0Val) + (vaVal) * std::sin (2.0*mpi*((std::real(phaseVal))/360)) ;
      }
      else
      {
        // 2PI to convert from hz to radians/sec
        SourceValue = (v0Val) + (vaVal) * std::sin (2.0*mpi*((std::real(freqVal))*std::real(time) + (std::real(phaseVal))/360)) * std::exp( -(std::real(time)*(std::real(thetaVal))));
      }
      return SourceValue;
    }

    // Note: this is only set up to compute dx w.r.t. time, for supporting breakpoints.  
    // And, it assumes that "time" is a special Op, and can't be a more complicated expression.
    // This function should be expanded to compute derivatives for other input params, such as va, v0, etc.
    virtual ScalarT dx (int i)
    {
      ScalarT dSource_dt = 0.0;

      ScalarT dTime_dt = this->time_->dx(i); 
      if(std::real(dTime_dt) != 0.0)
      { 
        ScalarT time = std::real(this->time_->val());
        ScalarT tdVal = this->childrenAstNodes_[3]->val();
        time -= std::real(tdVal);
        double mpi = M_PI;

        if (std::real(time) <= 0)
        {
          dSource_dt = 0.0;
        }
        else
        {
          // time derivative computed via Maple:
          std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
          //ScalarT v0Val    = child[0]->val();
          ScalarT vaVal    = child[1]->val();
          ScalarT freqVal  = child[2]->val();
          //ScalarT tdVal    = child[3]->val();
          ScalarT thetaVal = child[4]->val();
          ScalarT phaseVal = child[5]->val();

          dSource_dt = 2.0*vaVal*mpi*std::real(freqVal)*cos(2.0*mpi*(std::real(freqVal)*time+1/360*std::real(phaseVal)))*exp(-time*std::real(thetaVal))-vaVal*sin(2.0*mpi*(std::real(freqVal)*time+1/360*std::real(phaseVal)))*std::real(thetaVal)*exp(-time*std::real(thetaVal));
          dSource_dt *= dTime_dt;
        }
      }
      return dSource_dt;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) 
    {
      ScalarT dSource_dt = 0.0;
      if (!freqGiven_ && finalTime_ != 0.0)  
      {
        Teuchos::RCP<numval<ScalarT> > freqTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (this->childrenAstNodes_[2]);
        freqTmpOp->number = 1.0/finalTime_;
      }

      //ScalarT time = std::real(this->time_->val());
      ScalarT time;
      this->time_->dx2(time,derivs); // ERK check this!
      time = std::real(time);

      ScalarT tdVal = this->childrenAstNodes_[3]->val();
      time -= std::real(tdVal);
      double mpi = M_PI;

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      ScalarT v0Val    = child[0]->val();
      ScalarT vaVal    = child[1]->val();
      ScalarT freqVal  = child[2]->val();
      //ScalarT tdVal    = child[3]->val();
      ScalarT thetaVal = child[4]->val();
      ScalarT phaseVal = child[5]->val();

      if (std::real(time) <= 0)
      {
        dSource_dt = 0.0;
        result = (v0Val) + (vaVal) * std::sin (2.0*mpi*((std::real(phaseVal))/360)) ;
        if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
      }
      else
      {
        // 2PI to convert from hz to radians/sec
        result = (v0Val) + (vaVal) * std::sin (2.0*mpi*((std::real(freqVal))*std::real(time) + (std::real(phaseVal))/360)) * std::exp( -(std::real(time)*(std::real(thetaVal))));

          // time derivative computed via Maple:
        dSource_dt = 2.0*vaVal*mpi*std::real(freqVal)*cos(2.0*mpi*(std::real(freqVal)*time+1/360*std::real(phaseVal)))*exp(-time*std::real(thetaVal))-vaVal*sin(2.0*mpi*(std::real(freqVal)*time+1/360*std::real(phaseVal)))*std::real(thetaVal)*exp(-time*std::real(thetaVal));
      
        int numDerivs = derivs.size();
        for (int i=0;i<numDerivs;i++)
        {
          derivs[i] *= dSource_dt;
        }
      }
    }

    // in practice, only used for transient. update if this changes.
    virtual bool getIsComplex () { return false; }

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      if (tdGiven_)
      {
        double basetime=0.0;
        ScalarT tdVal = this->childrenAstNodes_[3]->val();
        ScalarT TD = std::real(tdVal);
        breakPointTimes.push_back(std::real(basetime+TD));
      }
      return true;
    }

    virtual void setFinalTime(double finalTime) { finalTime_ = finalTime; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice sin operator id = " << this->id_ << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[0]);
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[1]);
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[2]);
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[3]);
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[4]);
      AST_CALL_SUBOUTPUT(this->childrenAstNodes_[5]);
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "spice sin operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_sin codeGen function is not implemented yet" <<std::endl;
    }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool srcType() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      this->thisAstNode_ = thisAst_;
      Teuchos::RCP<spiceSinOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<spiceSinOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch
                                   
      this->childrenAstNodes_[0]->accept(visitor,this->childrenAstNodes_[0]);
      this->childrenAstNodes_[1]->accept(visitor,this->childrenAstNodes_[1]);
      this->childrenAstNodes_[2]->accept(visitor,this->childrenAstNodes_[2]);
      this->childrenAstNodes_[3]->accept(visitor,this->childrenAstNodes_[3]);
      this->childrenAstNodes_[4]->accept(visitor,this->childrenAstNodes_[4]);
      this->childrenAstNodes_[5]->accept(visitor,this->childrenAstNodes_[5]);

      time_->accept(visitor,time_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > time_;
    bool v0Given_, vaGiven_, freqGiven_, tdGiven_, thetaGiven_, phaseGiven_;
    double finalTime_;
};

//-------------------------------------------------------------------------------
// spice exp operator
//
// exp (v1,v2,td1,tau1,td2,tau2)
//
template <typename ScalarT>
class spiceExpOp : public astNode<ScalarT>
{
  public:
    spiceExpOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(args), 
      time_(time),
      v1Given_(false), v2Given_(false), td1Given_(false), tau1Given_(false), td2Given_(false), tau2Given_(false),
      startingTimeStep_(0.0)
    {
      if (args.size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_exp) needs at least 2 argument.  V1 and V2 are required for the EXP source function.")); yyerror(errStr);
      }

      if (args.size() > 6)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_exp) has too many arguments")); yyerror(errStr);
      }
   
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      if (child.size() < 6) { child.resize(6); }

      if (args.size() >= 1) { v1Given_=true;    } else { child[0] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 2) { v2Given_=true;    } else { child[1] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 3) { td1Given_=true;   } else { child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 4) { tau1Given_=true;  } else { child[3] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 5) { td2Given_=true;   } else { child[4] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 6) { tau2Given_=true;  } else { child[5] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    };

    virtual ScalarT val()
    {
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & v1_   = child[0];
      Teuchos::RCP<astNode<ScalarT> > & v2_   = child[1];
      Teuchos::RCP<astNode<ScalarT> > & td1_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & tau1_ = child[3];
      Teuchos::RCP<astNode<ScalarT> > & td2_  = child[4];
      Teuchos::RCP<astNode<ScalarT> > & tau2_ = child[5];

      // If neccessary, set defaults:
      // double tstep = solState_.startingTimeStep_;
      // if (!TD1given)  TD1 = 0.0;
      // if (!TAU1given) TAU1 = tstep;
      // if (!TD2given)  TD2 = TD1 + tstep;
      // if (!TAU2given) TAU2 = tstep;
      if (!tau1Given_)  { Teuchos::RCP<numval<ScalarT> > tau1TmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (tau1_); tau1TmpOp->number = startingTimeStep_; }

      if (!td2Given_)   
      { 
        Teuchos::RCP<numval<ScalarT> > td1TmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (td1_); 
        Teuchos::RCP<numval<ScalarT> > td2TmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (td2_); 
        td2TmpOp->number = td1TmpOp->number + startingTimeStep_;
      }

      if (!tau2Given_)  { Teuchos::RCP<numval<ScalarT> > tau2TmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (tau2_); tau2TmpOp->number = startingTimeStep_; }
      ScalarT time = std::real(time_->val());
      ScalarT SourceValue = 0.0;

      ScalarT TD1 = std::real(td1_->val()), TD2 = std::real(td2_->val());
      if (std::real(time) <= std::real(TD1))
      {
        SourceValue = v1_->val();
      }
      else if (std::real(time) <= std::real(TD2))
      {
        ScalarT V1 = v1_->val(), V2 = v2_->val();
        ScalarT TAU1 = std::real(tau1_->val());
        SourceValue = V1 + (V2-V1)*(1.0-std::exp(-(std::real(time)-std::real(TD1))/std::real(TAU1)));
      }
      else
      {
        ScalarT V1 = v1_->val(), V2 = v2_->val();
        ScalarT TAU1 = std::real(tau1_->val()), TAU2 = std::real(tau2_->val());
        SourceValue = V1 + (V2-V1)*(1.0-std::exp(-(std::real(time)-std::real(TD1))/std::real(TAU1))) +
                           (V1-V2)*(1.0-std::exp(-(std::real(time)-std::real(TD2))/std::real(TAU2))) ;
      }
      return SourceValue;
    }

    virtual ScalarT dx (int i)
    {
      return 0.0;
    }

    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) 
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, only used for transient. update if this changes.
    virtual bool getIsComplex () { return false; }

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & td1_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & td2_  = child[4];

      if (td1Given_)
      {
        double basetime=0.0;
        ScalarT TD1 = std::real(td2_->val());
        breakPointTimes.push_back(std::real(basetime+TD1));
      }

      if (td2Given_)
      {
        double basetime=0.0;
        ScalarT TD2 = std::real(td2_->val());
        breakPointTimes.push_back(std::real(basetime+TD2));
      }
      return true;
    }

    virtual void setStartingTimeStep(double timeStep) { startingTimeStep_ = timeStep; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice exp operator id = " << this->id_ << std::endl;
      ++indent;

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & v1_   = child[0];
      Teuchos::RCP<astNode<ScalarT> > & v2_   = child[1];
      Teuchos::RCP<astNode<ScalarT> > & td1_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & tau1_ = child[3];
      Teuchos::RCP<astNode<ScalarT> > & td2_  = child[4];
      Teuchos::RCP<astNode<ScalarT> > & tau2_ = child[5];

      AST_CALL_SUBOUTPUT(v1_)
      AST_CALL_SUBOUTPUT(v2_)
      AST_CALL_SUBOUTPUT(td1_)
      AST_CALL_SUBOUTPUT(tau1_)
      AST_CALL_SUBOUTPUT(td2_)
      AST_CALL_SUBOUTPUT(tau2_)
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "spice exp operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_exp codeGen function is not implemented yet" <<std::endl;
    }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool srcType() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      this->thisAstNode_ = thisAst_;
      Teuchos::RCP<spiceExpOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<spiceExpOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & v1_   = child[0];
      Teuchos::RCP<astNode<ScalarT> > & v2_   = child[1];
      Teuchos::RCP<astNode<ScalarT> > & td1_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & tau1_ = child[3];
      Teuchos::RCP<astNode<ScalarT> > & td2_  = child[4];
      Teuchos::RCP<astNode<ScalarT> > & tau2_ = child[5];

      v1_->accept(visitor,v1_);
      v2_->accept(visitor,v2_);
      td1_->accept(visitor,td1_);
      tau1_->accept(visitor,tau1_);
      td2_->accept(visitor,td2_);
      tau2_->accept(visitor,tau2_);
      time_->accept(visitor,time_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > time_;
    bool v1Given_, v2Given_, td1Given_, tau1Given_, td2Given_, tau2Given_;
    double startingTimeStep_;
};

//-------------------------------------------------------------------------------
// spice exp operator
//
// sffm (v0,va,fc,mdi,fs)
//
template <typename ScalarT>
class spiceSffmOp : public astNode<ScalarT>
{
  public:
    spiceSffmOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > & args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(args), 
      time_(time),
      v0Given_(false), vaGiven_(false), fcGiven_(false), mdiGiven_(false), fsGiven_(false),
      finalTime_(0.0)
    {
      if (args.size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sffm) needs at least 2 argument.  V0 and VA are required for the SFFM source function.")); yyerror(errStr);
      }

      if (args.size() > 5)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sffm) has too many arguments")); yyerror(errStr);
      }
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      if (child.size() < 5) { child.resize(5); }

      if (args.size() >= 1) { v0Given_=true;    } else { child[0] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 2) { vaGiven_=true;    } else { child[1] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 3) { fcGiven_=true;    } else { child[2] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 4) { mdiGiven_=true;   } else { child[3] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args.size() >= 5) { fsGiven_=true;    } else { child[4] = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    };

    virtual ScalarT val()
    {
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & v0_  = child[0];
      Teuchos::RCP<astNode<ScalarT> > & va_  = child[1];
      Teuchos::RCP<astNode<ScalarT> > & fc_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & mdi_ = child[3];
      Teuchos::RCP<astNode<ScalarT> > & fs_  = child[4];

      // If neccessary, set the defaults:
      //double tstop = solState_.finalTime_;
      //if (!FCgiven) FC = 1.0/tstop;
      //if (!FSgiven) FS = 1.0/tstop;
      //
      if (!fcGiven_ && finalTime_ != 0.0) { Teuchos::RCP<numval<ScalarT> > fcTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (fc_); fcTmpOp->number = 1.0/finalTime_; }
      if (!fsGiven_ && finalTime_ != 0.0) { Teuchos::RCP<numval<ScalarT> > fsTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (fs_); fsTmpOp->number = 1.0/finalTime_; }

      ScalarT time = std::real(time_->val());

      ScalarT V0 = v0_->val();
      ScalarT VA = va_->val();
      ScalarT FC = fc_->val();
      ScalarT MDI = mdi_->val();
      ScalarT FS = fs_->val();

      double mpi = M_PI;
      ScalarT SourceValue = V0 + VA * sin((2 * mpi * FC * std::real(time)) +
                  MDI * sin (2 * mpi * FS * std::real(time)));

      return SourceValue;
    }

    virtual ScalarT dx (int i)
    {
      return 0.0;
    }
 
    virtual void dx2(ScalarT & result, std::vector<ScalarT> & derivs) 
    {
      result = val();
      if ( !(derivs.empty() ) ) { std::fill(derivs.begin(),derivs.end(),0.0);  }
    }

    // in practice, only used for transient. update if this changes.
    virtual bool getIsComplex () { return false; }

    virtual void setFinalTime(double finalTime) { finalTime_ = finalTime; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice sffm operator id = " << this->id_ << std::endl;
      ++indent;

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & v0_  = child[0];
      Teuchos::RCP<astNode<ScalarT> > & va_  = child[1];
      Teuchos::RCP<astNode<ScalarT> > & fc_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & mdi_ = child[3];
      Teuchos::RCP<astNode<ScalarT> > & fs_  = child[4];

      AST_CALL_SUBOUTPUT(v0_)
      AST_CALL_SUBOUTPUT(va_)
      AST_CALL_SUBOUTPUT(fc_)
      AST_CALL_SUBOUTPUT(mdi_)
      AST_CALL_SUBOUTPUT(fs_)
    }

    virtual void compactOutput(std::ostream & os)
    {
      os << "spice sffm operator id = " << this->id_ << std::endl;
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_sffm codeGen function is not implemented yet" <<std::endl;
    }

    virtual bool getIsTreeConstant() { return false; }
    virtual bool srcType() { return true; }

    virtual void accept (nodeVisitor<ScalarT> & visitor, Teuchos::RCP<astNode<ScalarT> > & thisAst_) 
    { 
      this->thisAstNode_ = thisAst_;
      Teuchos::RCP<spiceSffmOp<ScalarT> > castToThis = Teuchos::rcp_static_cast<spiceSffmOp<ScalarT> > (thisAst_);
      visitor.visit( castToThis ); // 2nd dispatch

      std::vector<Teuchos::RCP<astNode<ScalarT> > > & child =  this->childrenAstNodes_;
      Teuchos::RCP<astNode<ScalarT> > & v0_  = child[0];
      Teuchos::RCP<astNode<ScalarT> > & va_  = child[1];
      Teuchos::RCP<astNode<ScalarT> > & fc_  = child[2];
      Teuchos::RCP<astNode<ScalarT> > & mdi_ = child[3];
      Teuchos::RCP<astNode<ScalarT> > & fs_  = child[4];

      v0_->accept(visitor,v0_);
      va_->accept(visitor,va_);
      fc_->accept(visitor,fc_);
      mdi_->accept(visitor,mdi_);
      fs_->accept(visitor,fs_);
      time_->accept(visitor,time_);
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > time_;
    bool v0Given_, vaGiven_, fcGiven_, mdiGiven_, fsGiven_;
    double finalTime_;
};

#endif
