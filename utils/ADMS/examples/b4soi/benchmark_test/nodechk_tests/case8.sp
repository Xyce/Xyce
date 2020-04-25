*NODECHK Warnings
*CASE 8: DGSEPBT with TNODEOUT=1, SHMOD=0 
*Warning: "s2"
******************************************

.option abstol=1e-6 reltol=1e-6 node post probe ingold

.hdl "../../code/bsimsoi.va"
.include "modelcard.nmos"

* --- Voltage Sources ---
vd d  0 dc=0.05
vg g  0 dc=0
vs s  0 dc=0
ve e  0 dc=0
 
* --- Transistor ---
x1 d g s e p b t nmos1 W=10u L=0.5u SOIMOD=2 TNODEOUT=1 SHMOD=0 

* --- DC Analysis ---
*Id-Vg Characteristics for NMOS 
.dc vg 0 0.05 0.05
.print dc i(x1.d) 
.end
