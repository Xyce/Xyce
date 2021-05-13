% gilbert.m4
%
% 
%
.PS                            # Pic input begins with .PS
cct_init                       # Set defaults

elen = 0.75                    # Variables are allowed; default units are inches
Origin: Here                   # Position names are capitalized
  rgbdraw(0,0,1,GndV2: ground)
  V2: source(up_ elen,V); llabel(,V_2,); rlabel(,6V,);
  Lv2: line right_ 0.5*elen; dot; llabel(,\texttt{\color{red}1},);

  # first transistor pair (Q3,Q4)
  Rq3b: resistor(right_ 1.5*elen from Lv2.end); llabel(,R_5,); rlabel(,100 \Omega,); 
  dot; llabel(,\texttt{\color{red}2},);
  line 0.25*elen from Rq3b.end
  rgbdraw(1,0,1,Q3: bi_tr(up_ elen,L,N,E) with .B at Here); rlabel(,Q_3,)
  LQ3e2Q4e1: line right 0.5*elen from Q3.E;
  dot;
  LQ3e2Q4e2: line right 0.5*elen from LQ3e2Q4e1.end;
  rgbdraw(1,0,1,Q4: bi_tr(up_ elen,R,N,E) with .E at Here); llabel(,Q_4,)

  # base Q4 to base Q5 connection line
  LQ4b2Q5b1: line right 0.5*elen from Q4.B; dot; clabel(,,\texttt{\color{red}6}); 
  Point1: LQ4b2Q5b1.end;
  LQ4b2Q5b2: line right 0.5*elen from LQ4b2Q5b1.end;
  VLO: source(down_ 0.75*elen from Point1,S); llabel(+,,-);rlabel(,V_{LO},);

  # second transistor pair (Q5,Q6)
  rgbdraw(1,0,1,Q5: bi_tr(up_ elen,L,N,E) with .B at LQ4b2Q5b2.end); rlabel(,Q_5,)
  LQ5e2Q6e1: line right 0.5*elen from Q5.E;
  dot;
  LQ5e2Q6e2: line right 0.5*elen from LQ5e2Q6e1.end;
  rgbdraw(1,0,1,Q6: bi_tr(up_ elen,R,N,E) with .E at Here); llabel(,Q_6,)
  Lq6b: line 0.25*elen from Q6.B; dot; clabel(,,\texttt{\color{red}2});
  Rq6b: resistor(right_ 1.5*elen from Lq6b.end); llabel(,R_6,); rlabel(,100 \Omega,); 
 
  # the big line (conecting the base resistor of Q3 and base resistor of Q6 to V2=6.0V)
  L12a: line down_ 3.75*elen from Lv2.end;
  L12c: line down_ 3.75*elen from Rq6b.end;
  L12b: line from L12a.end to L12c.end;

  # connecting emitter of Q3 to collector of Q1, and setting up Q1:
  Lq3e2q1c: line down_ 0.5*elen from LQ3e2Q4e1.end; dot
  rgbdraw(1,0,1,Q1: bi_tr(up_ elen,L,N,E) with .C at Here); rlabel(,Q_1,)

  Lq1b: line left_ 0.25*elen from Q1.B; dot; rlabel(,\texttt{\color{red}10},);
  Rq1b: resistor(left_ elen from Lq1b.end); rlabel(,R_3,);llabel(,1500 \Omega,);
  Lq1b2: line down 0.75*elen from Rq1b.end; dot; rlabel(,,\texttt{\color{red}16});
  Lq1b3: line down 1.25*elen from Lq1b2.end;
  Lv1: line left_ 0.5*elen from Lq1b2.end;
  V1: reversed(`source',down_ elen,V); llabel(,,1.8V);rlabel(,V_1,);
  rgbdraw(0,0,1,GndV1: ground)

  # connecting emitter of Q6 to collector of Q2 :
  Lq6e2q2c: line down_ 0.5*elen from LQ5e2Q6e1.end; dot
  rgbdraw(1,0,1,Q2: bi_tr(up_ elen,R,N,E) with .C at Here); llabel(,Q_2,)

  Lq2b: line right_ 0.25*elen from Q2.B; dot; llabel(,\texttt{\color{red}15},);
  Rq2b: resistor(right_ elen from Lq2b.end); llabel(,R_4,);rlabel(,1500 \Omega,);
  Lq2b2: line down_ 2*elen from Rq2b.end; 

  # crossover line(s) connecting base of Q3 (inside base resistor) to base of Q6, and also the the VLO source.
  L222a: line down_ 0.75*elen from Rq3b.end;
  Point2: L222a.end;
  L222b: crossover(right_ from L222a.end to VLO.end,L,Lq3e2q1c); dot; rlabel(,\texttt{\color{red}2},);
  L222d: line down_ 0.75*elen from Lq6b.end;
  L222c: crossover(left_ from L222d.end to VLO.end,R,Lq6e2q2c); 


  # connect node 10 to node 15, with a Vin sin source in-between.  (big line, 3rd from bottom)
  Lnode10: line down_ 1.75*elen from Lq1b.end;
  Lnode15: line down_ 0.75*elen from Lq2b.end;
  V5: source(down_ elen,S); llabel(+,,-);rlabel(,V_5,);
  L10to15: line from V5.end to Lnode10.end;

  # connecting Q2 base resistor to the Q1 base resistor.  (big line, 2nd from bottom)
  Lq1b2q2b: line from Lq2b2.end to Lq1b3.end;

  # connecting Q1 RE to Q2 RE:
  dot( at Q1.E ); llabel(\texttt{\color{red}11},,);
  Rq1e: resistor(right_ from Q1.E to (Point1,Q1.E )); llabel(,R_1,);rlabel(,10 \Omega,);  dot; llabel(,\texttt{\color{red}12},);
  Rq2e: resistor(right_ from (Point1,Q1.E) to Q2.E); llabel(,R_2,);rlabel(,10 \Omega,);   dot; llabel(\texttt{\color{red}13},,);

  I1: source(down_ 0.75*elen from (Point1,Q1.E),I); rlabel(,I_1,); llabel(,1.8mA,);
  rgbdraw(0,0,1,GndI1: ground)


  # outputs:
  # Lout1
  Lq3c: line up 0.75*elen from Q3.C; dot; llabel(\texttt{\color{red}3},,);
  Lq5c: line up 0.75*elen from Q5.C; dot; llabel(\texttt{\color{red}3},,);
  Lout1a: line from Lq3c.end to Lq5c.end;

  
  # Lout2
  Lq4c: line up 0.25*elen from Q4.C; dot; llabel(\texttt{\color{red}5},,);
  Lq6c: line up 0.25*elen from Q6.C; dot; llabel(\texttt{\color{red}5},,);
  Lout2a: crossover(right_ from Lq4c.end to Lq6c.end,L,Lq5c);



  # R8 and Vdd (V3):
  R8: resistor(up_ 0.75*elen from Lq3c.end); llabel(,R_8,); rlabel(,1500 \Omega,);
  dot; rlabel(,,\texttt{\color{red}17});
  Lv3: line left_ elen from R8.end;
  V3: source(down_ elen from Lv3.end); llabel(+,,-);rlabel(V_{3},,8V);
  rgbdraw(0,0,1,ground)

  # R7:
  Lq6c2: line  up_ 0.5*elen from Lq6c.end;
  R7: resistor(up_ 0.75*elen from Lq6c2.end);llabel(,R_7,); rlabel(,1500 \Omega,);
  Lq6c2q3c: line left from R7.end to R8.end;


  #final Lout1 crossover:
  Lout1b: crossover(right_ 2*elen from Lout1a.end,L,Lq6c2);
#llabel(,\texttt{\color{red}3},);

  #final Lout2:
  Lout2b: line right_ elen from Lout2a.end;

.PE                            # Pic input ends
