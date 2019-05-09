
# 
# compareLatexTables.py new_dir
#
# this compares a set of tables used in Xyce's Reference Guide
# in this directory to some new ones generated in new_dir
# Some filtering is done before diff is called to allow for 
# different formatting commands
#
# this list of tables is fixed and must be edited in this
# file to change/add new tables

# this is a list of tuples of files to compare.  Since 
# the files created by Xyce have different names than 
# the ones used in the documentation we need this
# data structure to correlate one file with another.
# 
# Note: to generate Xyce's documentation tables, use
# Xyce -doc and Xyce -doc_cat

import sys
import os
import string

if __name__ == '__main__':
  
  # sys.argv[0] is the python program's file name
  newDirName = sys.argv[1]

  tableList = [ ( "C_1_Device_Params.tex",      "C_1_Device_Instance_Params.tex" ),
      ( "C_1_Model_Params.tex",                 "C_1_Device_Model_Params.tex" ),
      ( "D_1,2_Device_Params.tex",              "D_1,2_Device_Instance_Params.tex" ),
      ( "D_1,2_Model_Params.tex",               "D_1,2_Device_Model_Params.tex" ),
      ( "Digital_1_Device_Params.tex",          "Digital_1_Device_Instance_Params.tex" ),
      ( "Digital_1_Model_Params.tex",           "Digital_1_Device_Model_Params.tex" ),
      ( "J_1,2_Device_Params.tex",              "J_1,2_Device_Instance_Params.tex" ),
      ( "J_1,2_Model_Params.tex",               "J_1,2_Device_Model_Params.tex" ),
      ( "M_1_Device_Instance_Params.tex",       "M_1_Device_Instance_Params.tex" ),
      ( "M_1_Device_Instance_cat_Params.tex",   "M_1_Device_Instance_cat_Params.tex" ),
      ( "M_1_Device_Model_Params.tex",          "M_1_Device_Model_Params.tex" ),
      ( "M_1_Device_Model_cat_Params.tex",      "M_1_Device_Model_cat_Params.tex" ),
      ( "M_2_Device_Instance_Params.tex",       "M_2_Device_Instance_Params.tex" ),
      ( "M_2_Device_Instance_cat_Params.tex",   "M_2_Device_Instance_cat_Params.tex" ),
      ( "M_2_Device_Model_Params.tex",          "M_2_Device_Model_Params.tex" ),
      ( "M_2_Device_Model_cat_Params.tex",      "M_2_Device_Model_cat_Params.tex" ),
      ( "M_3_Device_Instance_Params.tex",       "M_3_Device_Instance_Params.tex" ),
      ( "M_3_Device_Instance_cat_Params.tex",   "M_3_Device_Instance_cat_Params.tex" ),
      ( "M_3_Device_Model_Params.tex",          "M_3_Device_Model_Params.tex" ),
      ( "M_3_Device_Model_cat_Params.tex",      "M_3_Device_Model_cat_Params.tex" ),
      ( "M_6_Device_Instance_Params.tex",       "M_6_Device_Instance_Params.tex" ),
      ( "M_6_Device_Instance_cat_Params.tex",   "M_6_Device_Instance_cat_Params.tex" ),
      ( "M_6_Device_Model_Params.tex",          "M_6_Device_Model_Params.tex" ),
      ( "M_6_Device_Model_cat_Params.tex",      "M_6_Device_Model_cat_Params.tex" ),
      ( "M_9_Device_Instance_Params.tex",       "M_9_Device_Instance_Params.tex" ),
      ( "M_9_Device_Instance_cat_Params.tex",   "M_9_Device_Instance_cat_Params.tex" ),
      ( "M_9_Device_Model_Params.tex",          "M_9_Device_Model_Params.tex" ),
      ( "M_9_Device_Model_cat_Params.tex",      "M_9_Device_Model_cat_Params.tex" ),
      ( "M_10_Device_Instance_Params.tex",      "M_10_Device_Instance_Params.tex" ), 
      ( "M_10_Device_Instance_cat_Params.tex",  "M_10_Device_Instance_cat_Params.tex" ),
      ( "M_10_Device_Model_Params.tex",         "M_10_Device_Model_Params.tex" ),
      ( "M_10_Device_Model_cat_Params.tex",     "M_10_Device_Model_cat_Params.tex" ),
      ( "M_14_Device_Instance_cat_Params.tex",  "M_14_Device_Instance_cat_Params.tex" ),
      ( "M_14_Device_Model_cat_Params.tex",     "M_14_Device_Model_cat_Params.tex" ),
      ( "M_18_Device_Params.tex",               "M_18_Device_Instance_Params.tex" ),
      ( "M_18_Model_Params.tex",                "M_18_Device_Model_Params.tex" ),
      ( "R_1_Device_Params.tex",                "R_1_Device_Instance_Params.tex" ),
      ( "R_1_Model_Params.tex",                 "R_1_Device_Model_Params.tex" ),
      ( "R_2_Device_Instance_Params.tex",       "R_2_Device_Instance_Params.tex" ),
      ( "R_2_Device_Model_Params.tex",          "R_2_Device_Model_Params.tex" ),
      ( "Z_1_Device_Params.tex",                "Z_1_Device_Instance_Params.tex" ),
      ( "Z_1_Model_Params.tex",                 "Z_1_Device_Model_Params.tex" ) ]

  # unknows matches
  # GStbl.tex
  # ICStbl.tex
  # IVStbl.tex
  # R1DeviceParams.tex
  # R1ModelParams.tex

  for tablePair in tableList:
    oldTable = tablePair[0]
    newTable = os.path.join(newDirName, tablePair[1])
    newTableText = newTable + ".text"
    oldTableText = oldTable + ".text"
    if( os.path.isfile( oldTable ) ):
      command = "grep -e 'tt' " + oldTable + " > " + oldTableText 
      os.system(command)
    if( os.path.isfile( oldTable ) ):
      command = "grep -e 'tt' " + newTable + " > " + newTableText
      os.system(command)
    diffResultsFileName = oldTable + "_diffResults.text"
    diffCommand = "diff -w " + oldTableText + " " + newTableText + " > " + diffResultsFileName
    diffResult = os.system( diffCommand )
    if( diffResult == 0 ):
      # no diff found so delete this file
      # clean up temp files
      os.remove( diffResultsFileName )
      os.remove( oldTableText )
      os.remove( newTableText )
    else:
      print "Difference found between ", oldTable, " and ", newTable
      # we don't clean up old/new table text files in this case
     


