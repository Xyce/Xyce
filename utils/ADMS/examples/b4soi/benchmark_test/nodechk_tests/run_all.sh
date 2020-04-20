#!/bin/bash
hspice -hdlpath ../../Code/ -i case1.sp -o case1
hspice -hdlpath ../../Code/ -i case2.sp -o case2
hspice -hdlpath ../../Code/ -i case3.sp -o case3
hspice -hdlpath ../../Code/ -i case4.sp -o case4
hspice -hdlpath ../../Code/ -i case5.sp -o case5
hspice -hdlpath ../../Code/ -i case6.sp -o case6
hspice -hdlpath ../../Code/ -i case7.sp -o case7
hspice -hdlpath ../../Code/ -i case8.sp -o case8
hspice -hdlpath ../../Code/ -i case9.sp -o case9
hspice -hdlpath ../../Code/ -i case10.sp -o case10
hspice -hdlpath ../../Code/ -i case11.sp -o case11
hspice -hdlpath ../../Code/ -i case12.sp -o case12
hspice -hdlpath ../../Code/ -i case13.sp -o case13
hspice -hdlpath ../../Code/ -i case14.sp -o case14
hspice -hdlpath ../../Code/ -i case15.sp -o case15
hspice -hdlpath ../../Code/ -i case16.sp -o case16

rm -f *.ac* *.ic* *.pa* *.st* *.val *.sw* *.tr* *.sc* *.mt* *.ma* sxcmd.log acsymm.sp acsymm.lis *.out *.ms* *.valog
rm -rf *.pvadir *.ahdlSimDB
