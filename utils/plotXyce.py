#!/usr/bin/env python
from __future__ import print_function
#-------------------------------------------------------------------------------
#
# File: plotXyce.py
#
# Purpose: To quickly plot and/or animate Xyce output files in matlab style.
#
#          This plotting program suppots Xyce's 3 major text 
#          file formats: std, tecplot, and probe (csd).
#
#          This plotting program even supports .step output.
#
# Usage:  plotXyce.py [options] foo.cir.prn
#
# Author: Todd Coffey, SNL/NM
#
#-------------------------------------------------------------------------------
"""
This script plots standard Xyce output files
If there are 4 or fewer variables to plot, they will be put in one figure,
otherwise each variable will get its own figure.

Usage:  plotXyce.py [options] foo.cir.prn
options:
  -h or --help                    this display
  -v or --verbose                 print verbose output
  -p plotopts or --plot=plotopts  add plotoptions to plot commands
  -f or --figures                 put each veriable in its own figure

plotopts can be the following and can be repeated:
linestyles: - -- -. : . , o ^ v < > s + x D d 1 2 3 4 h H p | _ 
colors: b g r c m y k w
linewidth=2
"""
#Read:
#Python Essential Reference by Beazley

import os, sys, re
import numpy
from findBlock import findBlock

#-------------------------------------------------------------------------------
def getXyceData(file,verbose=False):
  """
  getdata(file) reads the data from a Xyce output prn file into a list of
  column names and an array containing the data.
  """
  if os.path.exists(file):
    input = open(file,'r').readlines()
    numlines = len(input)-2
    # remove spaces and braces from expressions
    line0 = input[0]
    match = findBlock(line0,delim="\{")
    while match[0]:
      beg = match[0]
      end = match[1]
      group = line0[beg+1:end-1] # remove braces
      line0 = line0[:beg] + re.sub(r"[ ]",r"",group) + line0[end+1:]
      match = findBlock(line0,end,delim="\{")
    tags = line0.split()
    numentries = len(tags)
    if verbose:
      print ("Read file: " + file)
      print ("Found columns:  ", end=' ')
      print (tags)
      print ("Found " + str(numlines) + " lines of data" )
    data = numpy.zeros((numlines,numentries),'double')
    if verbose:
      print ("Reading lines: ",)
    for i in range(1,len(input)-1):
      if verbose:
        print (".", end=' ')
      s = input[i].split()
      if s[0] != 'End':  # if this is the final line text string, skip it.
        if len(s) != numentries:
          print ("Error: " + file + ":" + str(i+1))
          print ("Number of columns read is not equal to number of columns in header.")
          print ("numentries = " + str(numentries))
          print ("len(s) = " + str(len(s)))
          sys.exit(1)
        data[i-1,:] = [float(j) for j in s]
    if verbose:
      print ("\n")
  else:
    print ("Error, file, " + file + " does not exist")
    print (getXyceData.__doc__)
    sys.exit(1)
  return (tags,data)

#-------------------------------------------------------------------------------
def parsePlotOpts(plotopts,verbose=False):
  """ Parse plot options into a list and a dictionary
  """
  if verbose and plotopts:
    print ("Plotting options = ",end=' ')
    print (plotopts)
  listopts = []
  dictopts = {}
  for j in range(len(plotopts)):
    if "=" in plotopts[j]:
      s = plotopts[j].split("=")
      dictopts[s[0]] = s[1]
    else:
      listopts.append(plotopts[j])
  if verbose:
    if listopts:
      print ("listopts = ", end=' ') ; print (listopts)
    if dictopts:
      print ("dictopts = ", end=' ') ; print (dictopts)
  return (listopts,dictopts)

#-------------------------------------------------------------------------------
import pylab
def plotXyceData(tags,data,file,verbose=False,figures=False,plotopts=[],showplot=True):
  """
  plotdata(tags,data,file,options) plots the data from a Xyce prn file labeling the
  plots with the column tags from the prn file and the name of the file that was
  read in.
  verbose     # verbosity [default=False]
  figures     # use separate figures for each plot [default=False]
  plotopts    # specific plot options [default=none]
  showplot    # call pylab.show() at end [default=True]
  """
  numentries = len(tags)
  dosubplot=True
  if "INDEX" in tags[0].upper(): # This is to handle xyce_verify.pl plotfile output
    indep=1
    numfigs=len(tags)-2
  else:
    indep=0
    numfigs=len(tags)-1
  if numentries > 6 or figures:
    dosubplot=False
  if verbose:
    print ("Number of columns to plot = " + str(numfigs))
    if dosubplot:
      print ("Using subplots for plotting")
    else:
      print ("Using individual figures for each variable")
  pylab.figure()
  for i in range(numfigs):
    if verbose:
      print ("Plotting column " + str(i+1) + " = " + tags[i+indep+1])
    if dosubplot:
      pylab.subplot(numfigs, 1, i+1)
    else:
      if i > 0: 
        pylab.figure()
    listopts,dictopts = parsePlotOpts(plotopts,verbose)
    pylab.plot(data[:,indep],data[:,i+indep+1],*listopts,**dictopts)
    pylab.ylabel(tags[i+indep+1])
    if i==0 or not dosubplot:
      pylab.title(file)
    if i==numfigs-1 or not dosubplot:
      pylab.xlabel(tags[indep])
  if showplot:
    pylab.show()


#-------------------------------------------------------------------------------
from getopt import getopt
def main():
  """
  This script plots standard Xyce output files
  If there are 4 or fewer variables to plot, they will be put in one figure,
  otherwise each variable will get its own figure.

  Usage:  plotXyce.py [options] foo.cir.prn
  options:
    -h or --help                    this display
    -v or --verbose                 print verbose output
    -p plotopts or --plot=plotopts  add plotoptions to plot commands
    -f or --figures                 put each veriable in its own figure

  plotopts can be the following and can be repeated:
  linestyles: - -- -. : . , o ^ v < > s + x D d 1 2 3 4 h H p | _ 
  colors: b g r c m y k w
  linewidth=2
  """
  poptions = { "plotopts":[], "verbose":False, "figures":0 }
  progDir,progName = os.path.split(sys.argv[0])
  options = "hvp:f"
  long_options = ["help","verbose","plot=","figures"]
  try:
    opts,args = getopt(sys.argv[1:],options,long_options)
  except:
    print ("Unrecognized argument")
    print (main.__doc__)
    sys.exit(1)
  for flag in opts:
    if flag[0] in ("-h","--help"):
      print (main.__doc__)
      sys.exit()
    elif flag[0] in ("-p","--plot"):
      poptions["plotopts"].append(flag[1])
    elif flag[0] in ("-v","--verbose"):
      poptions["verbose"] = True
    elif flag[0] in ("-f","--figures"):
      poptions["figures"] = 1
    else:
      print ("Unrecognized flag:", flag[0])
      print (main.__doc__)
      sys.exit(1)
  if len(args)==0:
    print ("No prn file specified")
    print (main.__doc__)
    sys.exit(1)
  tags,data = getXyceData(args[0],poptions["verbose"])
  plotXyceData(tags,data,args[0],**poptions)

if __name__ == "__main__":
  main()

