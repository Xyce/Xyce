*Sample netlist for BSIMSOI 
*Inverter DC

.option abstol=1e-6 reltol=1e-6 post

.hdl "../code/bsimsoi.va"
.include "modelcard.pmos"
.include "modelcard.nmos"

Vpower VD 0 1.5
Vgnd VS 0 0
Vgate Gate 0 0.0
xn0 VS Gate Out VS nmos1 W=10u L=0.18u
xp0 VD Gate Out VS pmos1 W=20u L=0.18u

.dc Vgate 0 1.5 0.05
.print dc v(out)
.END
