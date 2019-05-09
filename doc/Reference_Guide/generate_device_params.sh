#!/bin/sh

XYCE=~/devl/Xyce/build/Linux-x86_64-serial/src/Xyce

$XYCE -doc_cat M 10
$XYCE -doc_cat M 14
$XYCE -doc_cat M 9
$XYCE -doc_cat M 1
$XYCE -doc_cat M 2
$XYCE -doc_cat M 3
$XYCE -doc_cat M 6


$XYCE -doc C 1
$XYCE -doc D 1
$XYCE -doc D 2
$XYCE -doc J 1
$XYCE -doc J 2
$XYCE -doc L 1
$XYCE -doc M 103
$XYCE -doc M 18
$XYCE -doc M 1
$XYCE -doc M 2
$XYCE -doc M 301
$XYCE -doc M 3
$XYCE -doc M 6
$XYCE -doc K 1
$XYCE -doc O 1
$XYCE -doc Q 10
$XYCE -doc Q 1
$XYCE -doc Q 23
$XYCE -doc R 1
$XYCE -doc R 2
$XYCE -doc S 1
$XYCE -doc T 1
$XYCE -doc U 1
$XYCE -doc YNot 1
$XYCE -doc YMin 1
$XYCE -doc Z 1

