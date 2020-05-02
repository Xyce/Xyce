
#ifndef _ast_spice_src_h_
#define _ast_spice_src_h_

#define AST_CALL_SUBFUNC(PTR,FUNC,ARG)  if(this->PTR) { this->PTR->FUNC(ARG);  }

#define AST_CALL_SUBOUTPUT(PTR)  if(  !(Teuchos::is_null(this->PTR))  ) { os << std::setw(indent) << " "; os << #PTR << ": " << std::endl; this->PTR->output(os,indent+1); }

//-------------------------------------------------------------------------------
// spice pulse  operator
//
// pulse( v1,v2,td,tr,tf,pw,per)
//
template <typename ScalarT>
class spicePulseOp : public astNode<ScalarT>
{
  public:
    spicePulseOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > * args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(), 
      v1Given_(false), v2Given_(false), tdGiven_(false),
      trGiven_(false), tfGiven_(false), pwGiven_(false), perGiven_(false),
      time_(time), bpTol_(0.0), startingTimeStep_(0.0), finalTime_(0.0)
  {
    if (args->size() < 1)
    {
      std::vector<std::string> errStr(1,std::string("AST node (spice_pulse) needs at least 1 argument.  V1 is required for the PULSE source function.")); yyerror(errStr);
    }

    if (args->size() > 7)
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
    if (args->size() >= 1) { v1_ = (*args)[0]; v1Given_=true; } else { v1_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    if (args->size() >= 2) { v2_ = (*args)[1]; v2Given_=true; } else { v2_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    if (args->size() >= 3) { td_ = (*args)[2]; tdGiven_=true; } else { td_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    if (args->size() >= 4) { tr_ = (*args)[3]; trGiven_=true; } else { tr_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    if (args->size() >= 5) { tf_ = (*args)[4]; tfGiven_=true; } else { tf_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    if (args->size() >= 6) { pw_ = (*args)[5]; pwGiven_=true; } else { pw_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    if (args->size() >= 7) { per_ = (*args)[6]; perGiven_=true; } else { per_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
  };

    virtual ScalarT val()
    {
      if (!trGiven_) { Teuchos::RCP<numval<ScalarT> > trTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (tr_); trTmpOp->number = startingTimeStep_; }
      if (!tfGiven_) { Teuchos::RCP<numval<ScalarT> > tfTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (tf_); tfTmpOp->number = startingTimeStep_; }
      if (!pwGiven_) { Teuchos::RCP<numval<ScalarT> > pwTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (pw_); pwTmpOp->number = finalTime_; }
      if (!perGiven_) { Teuchos::RCP<numval<ScalarT> > perTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (per_); perTmpOp->number = finalTime_; }

      ScalarT time = std::real(this->time_->val());
      ScalarT V1 = this->v1_->val();
      ScalarT V2 = this->v2_->val();
      ScalarT TD = std::real(this->td_->val());
      ScalarT TR = std::real(this->tr_->val());
      ScalarT TF = std::real(this->tf_->val());
      ScalarT PW = std::real(this->pw_->val());
      ScalarT PER = std::real(this->per_->val());

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

    virtual bool getBreakPoints(std::vector<Xyce::Util::BreakPoint> & breakPointTimes)
    {
      ScalarT time = std::real(this->time_->val());
      ScalarT V1 = this->v1_->val();
      ScalarT V2 = this->v2_->val();
      ScalarT TD = std::real(this->td_->val());
      ScalarT TR = std::real(this->tr_->val());
      ScalarT TF = std::real(this->tf_->val());
      ScalarT PW = std::real(this->pw_->val());
      ScalarT PER = std::real(this->per_->val());

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
      os << "spice pulse operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v1_)
      AST_CALL_SUBOUTPUT(v2_)
      AST_CALL_SUBOUTPUT(td_)
      AST_CALL_SUBOUTPUT(tr_)
      AST_CALL_SUBOUTPUT(tf_)
      AST_CALL_SUBOUTPUT(pw_)
      AST_CALL_SUBOUTPUT(per_)
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_pulse codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(v1_) AST_GET_INTERESTING_OPS(v2_) AST_GET_INTERESTING_OPS(td_)
AST_GET_INTERESTING_OPS(tr_) AST_GET_INTERESTING_OPS(tf_) AST_GET_INTERESTING_OPS(pw_)
AST_GET_INTERESTING_OPS(per_) AST_GET_INTERESTING_OPS(time_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(v1_) AST_GET_PARAM_OPS(v2_) AST_GET_PARAM_OPS(td_)
AST_GET_PARAM_OPS(tr_) AST_GET_PARAM_OPS(tf_) AST_GET_PARAM_OPS(pw_)
AST_GET_PARAM_OPS(per_) AST_GET_PARAM_OPS(time_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(v1_) AST_GET_FUNC_ARG_OPS(v2_) AST_GET_FUNC_ARG_OPS(td_)
AST_GET_FUNC_ARG_OPS(tr_) AST_GET_FUNC_ARG_OPS(tf_) AST_GET_FUNC_ARG_OPS(pw_)
AST_GET_FUNC_ARG_OPS(per_) AST_GET_FUNC_ARG_OPS(time_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(v1_) AST_GET_FUNC_OPS(v2_) AST_GET_FUNC_OPS(td_)
AST_GET_FUNC_OPS(tr_) AST_GET_FUNC_OPS(tf_) AST_GET_FUNC_OPS(pw_)
AST_GET_FUNC_OPS(per_) AST_GET_FUNC_OPS(time_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(v1_) AST_GET_VOLT_OPS(v2_) AST_GET_VOLT_OPS(td_)
AST_GET_VOLT_OPS(tr_) AST_GET_VOLT_OPS(tf_) AST_GET_VOLT_OPS(pw_)
AST_GET_VOLT_OPS(per_) AST_GET_VOLT_OPS(time_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(v1_) AST_GET_CURRENT_OPS(v2_) AST_GET_CURRENT_OPS(td_)
AST_GET_CURRENT_OPS(tr_) AST_GET_CURRENT_OPS(tf_) AST_GET_CURRENT_OPS(pw_)
AST_GET_CURRENT_OPS(per_) AST_GET_CURRENT_OPS(time_)
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > v1_, v2_, td_, tr_, tf_, pw_, per_, time_;
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
    spiceSinOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > * args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(), 
      time_(time),
      v0Given_(false), vaGiven_(false), freqGiven_(false), tdGiven_(false), thetaGiven_(false), phaseGiven_(false),
      finalTime_(0.0)
    {
      if (args->size() < 3)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sin) needs at least 3 argument.  V0, VA and FREQ are required for the SIN source function.")); yyerror(errStr);
      }

      if (args->size() > 6)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sin) has too many arguments")); yyerror(errStr);
      } 
    
      if (args->size() >= 1) { v0_    = (*args)[0]; v0Given_=true;    } else { v0_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 2) { va_    = (*args)[1]; vaGiven_=true;    } else { va_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 3) { freq_  = (*args)[2]; freqGiven_=true;  } else { freq_  = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 4) { td_    = (*args)[3]; tdGiven_=true;    } else { td_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 5) { theta_ = (*args)[4]; thetaGiven_=true; } else { theta_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 6) { phase_ = (*args)[5]; phaseGiven_=true; } else { phase_ = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    };

    virtual ScalarT val()
    {
      // defaults are as follows: (have to be set here b/c in constructor not sure if we know the final time yet
      //
      //   double tstop = solState_.finalTime_;
      //  if (!FREQgiven)  FREQ  = 1.0/tstop;

      if (!freqGiven_ && finalTime_ != 0.0)  
      {
        Teuchos::RCP<numval<ScalarT> > freqTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (freq_);
        freqTmpOp->number = 1.0/finalTime_;
      }

      ScalarT time = std::real(this->time_->val());
      time -= std::real(this->td_->val());
      double mpi = M_PI;
      ScalarT SourceValue = 0.0;

      if (std::real(time) <= 0)
      {
        SourceValue = (this->v0_->val()) + (this->va_->val()) * std::sin (2.0*mpi*((std::real(this->phase_->val()))/360)) ;
      }
      else
      {
        // 2PI to convert from hz to radians/sec
        SourceValue = (this->v0_->val()) + (this->va_->val()) * std::sin (2.0*mpi*((std::real(this->freq_->val()))*std::real(time) + (std::real(this->phase_->val()))/360)) * std::exp( -(std::real(time)*(std::real(this->theta_->val()))));
      }
      return SourceValue;
    }

    virtual ScalarT dx (int i)
    {
      return 0.0;
    }

    virtual void setFinalTime(double finalTime) { finalTime_ = finalTime; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice sin operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v0_)
      AST_CALL_SUBOUTPUT(va_)
      AST_CALL_SUBOUTPUT(freq_)
      AST_CALL_SUBOUTPUT(td_)
      AST_CALL_SUBOUTPUT(theta_)
      AST_CALL_SUBOUTPUT(phase_)
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_sin codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(v0_) AST_GET_INTERESTING_OPS(va_) AST_GET_INTERESTING_OPS(freq_)
AST_GET_INTERESTING_OPS(td_) AST_GET_INTERESTING_OPS(theta_) AST_GET_INTERESTING_OPS(phase_) AST_GET_INTERESTING_OPS(time_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(v0_) AST_GET_PARAM_OPS(va_) AST_GET_PARAM_OPS(freq_)
AST_GET_PARAM_OPS(td_) AST_GET_PARAM_OPS(theta_) AST_GET_PARAM_OPS(phase_) AST_GET_PARAM_OPS(time_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(v0_) AST_GET_FUNC_ARG_OPS(va_) AST_GET_FUNC_ARG_OPS(freq_)
AST_GET_FUNC_ARG_OPS(td_) AST_GET_FUNC_ARG_OPS(theta_) AST_GET_FUNC_ARG_OPS(phase_) AST_GET_FUNC_ARG_OPS(time_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(v0_) AST_GET_FUNC_OPS(va_) AST_GET_FUNC_OPS(freq_)
AST_GET_FUNC_OPS(td_) AST_GET_FUNC_OPS(theta_) AST_GET_FUNC_OPS(phase_) AST_GET_FUNC_OPS(time_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(v0_) AST_GET_VOLT_OPS(va_) AST_GET_VOLT_OPS(freq_)
AST_GET_VOLT_OPS(td_) AST_GET_VOLT_OPS(theta_) AST_GET_VOLT_OPS(phase_) AST_GET_VOLT_OPS(time_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(v0_) AST_GET_CURRENT_OPS(va_) AST_GET_CURRENT_OPS(freq_)
AST_GET_CURRENT_OPS(td_) AST_GET_CURRENT_OPS(theta_) AST_GET_CURRENT_OPS(phase_) AST_GET_CURRENT_OPS(time_)
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > v0_, va_, freq_, td_, theta_, phase_, time_;
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
    spiceExpOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > * args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(), 
      time_(time),
      v1Given_(false), v2Given_(false), td1Given_(false), tau1Given_(false), td2Given_(false), tau2Given_(false),
      startingTimeStep_(0.0)
    {
      if (args->size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_exp) needs at least 2 argument.  V1 and V2 are required for the EXP source function.")); yyerror(errStr);
      }

      if (args->size() > 6)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_exp) has too many arguments")); yyerror(errStr);
      } 
    
      if (args->size() >= 1) { v1_    = (*args)[0]; v1Given_=true;    } else { v1_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 2) { v2_    = (*args)[1]; v2Given_=true;    } else { v2_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 3) { td1_   = (*args)[2]; td1Given_=true;   } else { td1_   = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 4) { tau1_  = (*args)[3]; tau1Given_=true;  } else { tau1_  = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 5) { td2_   = (*args)[4]; td2Given_=true;   } else { td2_   = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 6) { tau2_  = (*args)[5]; tau2Given_=true;  } else { tau2_  = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); } 
    };

    virtual ScalarT val()
    {
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

      ScalarT time = std::real(this->time_->val());
      ScalarT SourceValue = 0.0;

      ScalarT TD1 = std::real(this->td1_->val()), TD2 = std::real(this->td2_->val());
      if (std::real(time) <= std::real(TD1))
      {
        SourceValue = this->v1_->val();
      }
      else if (std::real(time) <= std::real(TD2))
      {
        ScalarT V1 = this->v1_->val(), V2 = this->v2_->val();
        ScalarT TAU1 = std::real(this->tau1_->val());
        SourceValue = V1 + (V2-V1)*(1.0-std::exp(-(std::real(time)-std::real(TD1))/std::real(TAU1)));
      }
      else
      {
        ScalarT V1 = this->v1_->val(), V2 = this->v2_->val();
        ScalarT TAU1 = std::real(this->tau1_->val()), TAU2 = std::real(this->tau2_->val());
        SourceValue = V1 + (V2-V1)*(1.0-std::exp(-(std::real(time)-std::real(TD1))/std::real(TAU1))) +
                           (V1-V2)*(1.0-std::exp(-(std::real(time)-std::real(TD2))/std::real(TAU2))) ;
      }
      return SourceValue;
    }

    virtual ScalarT dx (int i)
    {
      return 0.0;
    }

    virtual void setStartingTimeStep(double timeStep) { startingTimeStep_ = timeStep; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice exp operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v1_)
      AST_CALL_SUBOUTPUT(v2_)
      AST_CALL_SUBOUTPUT(td1_)
      AST_CALL_SUBOUTPUT(tau1_)
      AST_CALL_SUBOUTPUT(td2_)
      AST_CALL_SUBOUTPUT(tau2_)
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_exp codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(v1_) AST_GET_INTERESTING_OPS(v2_) AST_GET_INTERESTING_OPS(td1_)
AST_GET_INTERESTING_OPS(tau1_) AST_GET_INTERESTING_OPS(td2_) AST_GET_INTERESTING_OPS(tau2_) AST_GET_INTERESTING_OPS(time_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(v1_) AST_GET_PARAM_OPS(v2_) AST_GET_PARAM_OPS(td1_)
AST_GET_PARAM_OPS(tau1_) AST_GET_PARAM_OPS(td2_) AST_GET_PARAM_OPS(tau2_) AST_GET_PARAM_OPS(time_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(v1_) AST_GET_FUNC_ARG_OPS(v2_) AST_GET_FUNC_ARG_OPS(td1_)
AST_GET_FUNC_ARG_OPS(tau1_) AST_GET_FUNC_ARG_OPS(td2_) AST_GET_FUNC_ARG_OPS(tau2_) AST_GET_FUNC_ARG_OPS(time_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(v1_) AST_GET_FUNC_OPS(v2_) AST_GET_FUNC_OPS(td1_)
AST_GET_FUNC_OPS(tau1_) AST_GET_FUNC_OPS(td2_) AST_GET_FUNC_OPS(tau2_) AST_GET_FUNC_OPS(time_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(v1_) AST_GET_VOLT_OPS(v2_) AST_GET_VOLT_OPS(td1_)
AST_GET_VOLT_OPS(tau1_) AST_GET_VOLT_OPS(td2_) AST_GET_VOLT_OPS(tau2_) AST_GET_VOLT_OPS(time_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(v1_) AST_GET_CURRENT_OPS(v2_) AST_GET_CURRENT_OPS(td1_)
AST_GET_CURRENT_OPS(tau1_) AST_GET_CURRENT_OPS(td2_) AST_GET_CURRENT_OPS(tau2_) AST_GET_CURRENT_OPS(time_) 
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > v1_, v2_, td1_, tau1_, td2_, tau2_, time_;
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
    spiceSffmOp (std::vector<Teuchos::RCP<astNode<ScalarT> > > * args, Teuchos::RCP<astNode<ScalarT> > &time):
      astNode<ScalarT>(), 
      time_(time),
      v0Given_(false), vaGiven_(false), fcGiven_(false), mdiGiven_(false), fsGiven_(false),
      finalTime_(0.0)
    {
      if (args->size() < 2)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sffm) needs at least 2 argument.  V0 and VA are required for the SFFM source function.")); yyerror(errStr);
      }

      if (args->size() > 5)
      {
        std::vector<std::string> errStr(1,std::string("AST node (spice_sffm) has too many arguments")); yyerror(errStr);
      } 
    
      if (args->size() >= 1) { v0_    = (*args)[0]; v0Given_=true;    } else { v0_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 2) { va_    = (*args)[1]; vaGiven_=true;    } else { va_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 3) { fc_    = (*args)[2]; fcGiven_=true;    } else { fc_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 4) { mdi_   = (*args)[3]; mdiGiven_=true;   } else { mdi_   = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
      if (args->size() >= 5) { fs_    = (*args)[4]; fsGiven_=true;    } else { fs_    = Teuchos::RCP<astNode<ScalarT> >(new numval<ScalarT>(0.0)); }
    };

    virtual ScalarT val()
    {
      // If neccessary, set the defaults:
      //double tstop = solState_.finalTime_;
      //if (!FCgiven) FC = 1.0/tstop;
      //if (!FSgiven) FS = 1.0/tstop;
      //
      if (!fcGiven_ && finalTime_ != 0.0) { Teuchos::RCP<numval<ScalarT> > fcTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (fc_); fcTmpOp->number = 1.0/finalTime_; }
      if (!fsGiven_ && finalTime_ != 0.0) { Teuchos::RCP<numval<ScalarT> > fsTmpOp = Teuchos::rcp_static_cast<numval<ScalarT> > (fs_); fsTmpOp->number = 1.0/finalTime_; }

      ScalarT time = std::real(this->time_->val());

      ScalarT V0 = this->v0_->val();
      ScalarT VA = this->va_->val();
      ScalarT FC = this->fc_->val();
      ScalarT MDI = this->mdi_->val();
      ScalarT FS = this->fs_->val();

      double mpi = M_PI;
      ScalarT SourceValue = V0 + VA * sin((2 * mpi * FC * std::real(time)) +
                  MDI * sin (2 * mpi * FS * std::real(time)));

      return SourceValue;
    }

    virtual ScalarT dx (int i)
    {
      return 0.0;
    }

    virtual void setFinalTime(double finalTime) { finalTime_ = finalTime; }

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice sffm operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v0_)
      AST_CALL_SUBOUTPUT(va_)
      AST_CALL_SUBOUTPUT(fc_)
      AST_CALL_SUBOUTPUT(mdi_)
      AST_CALL_SUBOUTPUT(fs_)
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_sffm codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(opVectorContainers<ScalarT> & ovc)
    {
AST_GET_INTERESTING_OPS(v0_) AST_GET_INTERESTING_OPS(va_) AST_GET_INTERESTING_OPS(fc_)
AST_GET_INTERESTING_OPS(mdi_) AST_GET_INTERESTING_OPS(fs_) AST_GET_INTERESTING_OPS(time_)
    }

    virtual void getParamOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector)
    {
AST_GET_PARAM_OPS(v0_) AST_GET_PARAM_OPS(va_) AST_GET_PARAM_OPS(fc_)
AST_GET_PARAM_OPS(mdi_) AST_GET_PARAM_OPS(fs_) AST_GET_PARAM_OPS(time_)
    }

    virtual void getFuncArgOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcArgOpVector)
    {
AST_GET_FUNC_ARG_OPS(v0_) AST_GET_FUNC_ARG_OPS(va_) AST_GET_FUNC_ARG_OPS(fc_)
AST_GET_FUNC_ARG_OPS(mdi_) AST_GET_FUNC_ARG_OPS(fs_) AST_GET_FUNC_ARG_OPS(time_)
    }

    virtual void getFuncOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector)
    {
AST_GET_FUNC_OPS(v0_) AST_GET_FUNC_OPS(va_) AST_GET_FUNC_OPS(fc_)
AST_GET_FUNC_OPS(mdi_) AST_GET_FUNC_OPS(fs_) AST_GET_FUNC_OPS(time_)
    }

    virtual void getVoltageOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector)
    {
AST_GET_VOLT_OPS(v0_) AST_GET_VOLT_OPS(va_) AST_GET_VOLT_OPS(fc_)
AST_GET_VOLT_OPS(mdi_) AST_GET_VOLT_OPS(fs_) AST_GET_VOLT_OPS(time_)
    }

    virtual void getCurrentOps(std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
    {
AST_GET_CURRENT_OPS(v0_) AST_GET_CURRENT_OPS(va_) AST_GET_CURRENT_OPS(fc_)
AST_GET_CURRENT_OPS(mdi_) AST_GET_CURRENT_OPS(fs_) AST_GET_CURRENT_OPS(time_)
    }

  private:
    Teuchos::RCP<astNode<ScalarT> > v0_, va_, fc_, mdi_, fs_, time_;
    bool v0Given_, vaGiven_, fcGiven_, mdiGiven_, fsGiven_;
    double finalTime_;
};

#endif

