#!/bin/bash
hspice -hdlpath ../code -i test1.sp -o test1
hspice -hdlpath ../code -i test2.sp -o test2
hspice -hdlpath ../code -i test3.sp -o test3
hspice -hdlpath ../code -i test4.sp -o test4
hspice -hdlpath ../code -i test5.sp -o test5
hspice -hdlpath ../code -i test6.sp -o test6
hspice -hdlpath ../code -i test7.sp -o test7
hspice -hdlpath ../code -i test8.sp -o test8
hspice -hdlpath ../code -i inverter_dc.sp -o inverter_dc
hspice -hdlpath ../code -i inverter_tr.sp -o inverter_tr
hspice -hdlpath ../code -i ring_osc.sp -o ring_osc
      