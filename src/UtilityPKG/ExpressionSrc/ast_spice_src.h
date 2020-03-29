
#ifndef _ast_spice_src_h_
#define _ast_spice_src_h_

#define AST_CALL_SUBFUNC(PTR,FUNC,ARG)  if(this->PTR) { this->PTR->FUNC(ARG);  }

#define AST_CALL_SUBOUTPUT(PTR,OUTSTR)  if(  !(Teuchos::is_null(this->PTR))  ) { os << std::setw(indent) << " "; os << OUTSTR << ": " << std::endl; this->PTR->output(os,indent+1); }

#define AST_GET_INTERESTING_OPS(PTR) if( !(Teuchos::is_null(this->PTR)) ) {  \
  if (this->PTR->paramType()) { paramOpVector.push_back(this->PTR); }  \
  if (this->PTR->funcType())    { funcOpVector.push_back(this->PTR); } \
  if (this->PTR->voltageType()) { voltOpVector.push_back(this->PTR); } \
  if (this->PTR->currentType()) { currentOpVector.push_back(this->PTR); } \
  this->PTR->getInterestingOps(paramOpVector,funcOpVector,voltOpVector,currentOpVector); } 

#define AST_GET_PARAM_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->paramType()) { paramOpVector.push_back(this->PTR); } this->PTR->getParamOps(paramOpVector); }

#define AST_GET_FUNC_ARG_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->getFunctionArgType()) { funcArgOpVector.push_back(this->PTR); } this->PTR->getFuncArgOps(funcArgOpVector); }

#define AST_GET_FUNC_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->funcType()) { funcOpVector.push_back(this->PTR); } this->PTR->getFuncOps(funcOpVector); }

#define AST_GET_VOLT_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->voltageType()) { voltOpVector.push_back(this->PTR); } this->PTR->getVoltageOps(voltOpVector); }

#define AST_GET_CURRENT_OPS(PTR)  if( !(Teuchos::is_null(this->PTR)) ) { if (this->PTR->currentType()) { currentOpVector.push_back(this->PTR); } this->PTR->getCurrentOps(currentOpVector); }

//-------------------------------------------------------------------------------
// spice pulse  operator
//
// pulse( v1,v2,td,tr,tf,pw,per)
//
template <typename ScalarT>
class spicePulseOp : public astNode<ScalarT>
{
  public:
    spicePulseOp (
        Teuchos::RCP<astNode<ScalarT> > &v1, Teuchos::RCP<astNode<ScalarT> > &v2, Teuchos::RCP<astNode<ScalarT> > &td,
        Teuchos::RCP<astNode<ScalarT> > &tr, Teuchos::RCP<astNode<ScalarT> > &tf, Teuchos::RCP<astNode<ScalarT> > &pw,
        Teuchos::RCP<astNode<ScalarT> > &per, Teuchos::RCP<astNode<ScalarT> > &time
        ):
      astNode<ScalarT>(), v1_(v1), v2_(v2), td_(td), tr_(tr), tf_(tf), pw_(pw), per_(per), time_(time)
  {};

    virtual ScalarT val()
    {
      // If neccessary, set the defaults:
      //double tstep = solState_.startingTimeStep_;
      //double tstop = solState_.finalTime_;
      //if (!TDgiven)  TD  = 0.0;
      //if (!TRgiven)  TR  = tstep;
      //if (!TFgiven)  TF  = tstep;
      //if (!PWgiven)  PW  = tstop;
      //if (!PERgiven) PER = tstop;
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

      ScalarT bpTol = 0.0; // this needs to come from the group
      ScalarT SourceValue = 0.0;

      if (std::real(time) <= 0 || (std::real(time) > (std::real(TR) + std::real(PW) + std::real(TF)) &&
            (fabs (std::real(time) - std::real(TR+PW+TF)) > std::real(bpTol)) ) )
      {
        SourceValue = V1;
      }
      else if ((std::real(time) > std::real(TR) && fabs(std::real(time)-std::real(TR)) > std::real(bpTol))
          && (std::real(time) < (std::real(TR) + std::real(PW)) || fabs (std::real(time)-std::real(TR+PW))<std::real(bpTol)) )
      {
        SourceValue = V2;
      }
      else if (std::real(time) > 0 && (std::real(time) < std::real(TR) || fabs(std::real(time)-std::real(TR)) < std::real(bpTol)))
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

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice pulse operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v1_,"v1")
      AST_CALL_SUBOUTPUT(v2_,"v2")
      AST_CALL_SUBOUTPUT(td_,"td")
      AST_CALL_SUBOUTPUT(tr_,"tr")
      AST_CALL_SUBOUTPUT(tf_,"tf")
      AST_CALL_SUBOUTPUT(pw_,"pw")
      AST_CALL_SUBOUTPUT(per_,"per")
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_pulse codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
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
    spiceSinOp (
        Teuchos::RCP<astNode<ScalarT> > &v0, Teuchos::RCP<astNode<ScalarT> > &va, Teuchos::RCP<astNode<ScalarT> > &freq,
        Teuchos::RCP<astNode<ScalarT> > &td, Teuchos::RCP<astNode<ScalarT> > &theta, Teuchos::RCP<astNode<ScalarT> > &phase,
        Teuchos::RCP<astNode<ScalarT> > &time
        ):
      astNode<ScalarT>(), v0_(v0), va_(va), freq_(freq), td_(td), theta_(theta), phase_(phase), time_(time)
  {};

    virtual ScalarT val()
    {
      // ERK. do this somehow:
      //   double tstop = solState_.finalTime_;
      //  if (!FREQgiven)  FREQ  = 1.0/tstop;

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

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice sin operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v0_,"v0")
      AST_CALL_SUBOUTPUT(va_,"va")
      AST_CALL_SUBOUTPUT(freq_,"freq")
      AST_CALL_SUBOUTPUT(td_,"td")
      AST_CALL_SUBOUTPUT(theta_,"theta")
      AST_CALL_SUBOUTPUT(phase_,"phase")
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_sin codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
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
    spiceExpOp (
        Teuchos::RCP<astNode<ScalarT> > &v1, Teuchos::RCP<astNode<ScalarT> > &v2, Teuchos::RCP<astNode<ScalarT> > &td1,
        Teuchos::RCP<astNode<ScalarT> > &tau1, Teuchos::RCP<astNode<ScalarT> > &td2, Teuchos::RCP<astNode<ScalarT> > &tau2,
        Teuchos::RCP<astNode<ScalarT> > &time
        ):
      astNode<ScalarT>(), v1_(v1), v2_(v2), td1_(td1), tau1_(tau1), td2_(td2), tau2_(tau2), time_(time) {};

    virtual ScalarT val()
    {
      // If neccessary, set defaults:
      // double tstep = solState_.startingTimeStep_;
      // if (!TD1given)  TD1 = 0.0;
      // if (!TAU1given) TAU1 = tstep;
      // if (!TD2given)  TD2 = TD1 + tstep;
      // if (!TAU2given) TAU2 = tstep;

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

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice exp operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v1_,"v1")
      AST_CALL_SUBOUTPUT(v2_,"v2")
      AST_CALL_SUBOUTPUT(td1_,"td1")
      AST_CALL_SUBOUTPUT(tau1_,"tau1")
      AST_CALL_SUBOUTPUT(td2_,"td2")
      AST_CALL_SUBOUTPUT(tau2_,"tau2")
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_exp codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
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
AST_GET_CURRENT_OPS(tau1_) AST_GET_CURRENT_OPS(td2_) AST_GET_CURRENT_OPS(tau2_) AST_GET_CURRENT_OPS(time_) }

  private:
    Teuchos::RCP<astNode<ScalarT> > v1_, v2_, td1_, tau1_, td2_, tau2_, time_;
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
    spiceSffmOp (
        Teuchos::RCP<astNode<ScalarT> > &v0,
        Teuchos::RCP<astNode<ScalarT> > &va,
        Teuchos::RCP<astNode<ScalarT> > &fc,
        Teuchos::RCP<astNode<ScalarT> > &mdi,
        Teuchos::RCP<astNode<ScalarT> > &fs,
        Teuchos::RCP<astNode<ScalarT> > &time
        ):
      astNode<ScalarT>(), v0_(v0), va_(va), fc_(fc), mdi_(mdi), fs_(fs), time_(time) {};

    virtual ScalarT val()
    {
      // If neccessary, set the defaults:
      //double tstop = solState_.finalTime_;
      //if (!FCgiven) FC = 1.0/tstop;
      //if (!FSgiven) FS = 1.0/tstop;

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

    virtual void output(std::ostream & os, int indent=0)
    {
      os << std::setw(indent) << " ";
      os << "spice sffm operator " << std::endl;
      ++indent;

      AST_CALL_SUBOUTPUT(v0_,"v0")
      AST_CALL_SUBOUTPUT(va_,"va")
      AST_CALL_SUBOUTPUT(fc_,"fc")
      AST_CALL_SUBOUTPUT(mdi_,"mdi")
      AST_CALL_SUBOUTPUT(fs_,"fs")
    }

    virtual void codeGen (std::ostream & os )
    {
      os << "// spice_sffm codeGen function is not implemented yet" <<std::endl;
    }

    virtual void getInterestingOps(
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & paramOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & funcOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & voltOpVector,
      std::vector<Teuchos::RCP<astNode<ScalarT> > > & currentOpVector)
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
};

#endif

