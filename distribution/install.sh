#!/bin/sh

#   install.sh
#   Extracts tarball into user selected directory and configures shell script for
#   running Xyce from command line.  


echo "-----------------------------------------------------------------------------"
echo " Copyright Notice"
echo 
echo " Copyright (C) 2002-2024, Sandia Corporation, Albuquerque, NM, USA.  "
echo " See the output of Xyce -license for details. "
echo "-----------------------------------------------------------------------------"
echo
echo
echo "Preparing to install Xyce..."


#   setup path to needed progs
PATH="$PATH:/bin:/usr/bin:/usr/local/bin:/sbin:/usr/sbin"


# set envvars 
VER=5.3
#ORIG_PWD=`pwd`

ORIG_PWD=`dirname $0`
echo "ORIG_PWD is $ORIG_PWD"
if [ "x$ORIG_PWD" = "x." ]
then
   ORIG_PWD=`pwd`
fi

# trap errors and ctrl-c 
alert_user ()
{
  echo
  echo "Setup interrupted..."
  $CNOK && echo "Installation aborted:  alert_user() called" >> "$CONFIG_NOTES"
  exit 11;
}
trap "alert_user" 1 2 3 4 6 8 10 12 13 15


#   make sure we have the install tarball
if [ ! -f "$ORIG_PWD/xtl.tar.gz" ]; then
    echo "ERROR:  There are missing install files.  Please unpack a"
    echo "ERROR:  valid install tarball and run install.sh again."
    echo
    echo "ERROR:  You must run this install script from the"
    echo "ERROR:  directory where it is unpacked."
    echo
    echo "ERROR:  Installation aborted."
    exit 1;
fi


#   record configuration options
CONFIG_NOTES="$ORIG_PWD/install.log.$$"
echo "Xyce Installation Notes" > "$CONFIG_NOTES"


#   what plat is this?  (in case we have a mega install bundle :)
OS_NAME=`uname -s 2> /dev/null | tr "[:upper:]" "[:lower:]" 2> /dev/null`
echo "platform:  $OS_NAME" >> "$CONFIG_NOTES"


# set default prefix
INSTALLBASE="$HOME/Xyce" 
DONOTPROMPT=0

if [ $# = 0 ]
then
  echo "No command line args given"
else
  INSTALLBASE=$1
  DONOTPROMPT=1
fi
    
#   determine where to install Xyce    
INVALID="INVALID"
while [ "$INVALID" = "INVALID" ]
do
    if [ $DONOTPROMPT = 0 ]
    then
      # a ? means the last entry was bogus
      echo "Where should Xyce be installed?  [ $INSTALLBASE ]  "
      read ANS  
    
      #   user typed something so expand and check for validity
      if [ "x$ANS" != "x" ]; then
          INSTALLBASE=`eval echo $ANS`
      fi
    fi
   
    #   throw out bogus entries
    WALKER="DONE"
    case "$INSTALLBASE" in
        *\?* | *\** | *\$* )
            echo "WARNING:  $INSTALLBASE is not a valid path name."
	    echo
        ;;
        
	*\ * )
	    echo "WARNING:  Path must not contain spaces."
	;;

        * )
            WALKER="$INSTALLBASE"
        ;;
    esac
  
    #   make sure path is writeable

    while [ "$WALKER" != "DONE" ]
    do
        if [ -d "$WALKER" ]; then
            if [ ! -w "$WALKER" ]; then
                #   cannot write to this path so alert user and try again
                echo "WARNING:  $INSTALLBASE is not writeable."
                echo "WARNING:  Please select a new location."
                INSTALLBASE="?"
            else
                #   path is writeable and we can proceed to copy files
                INVALID="ok"
            fi
            #   reached root (/) so stop climbing
            WALKER="DONE" 
        else
            #   did not find an existing file so step up the hierarchy
            #   removing the trailing dirname and preceding slash
            WALKER=`echo $WALKER | sed 's/[/]$//' | sed 's/[^/]*$//'`
            
            #   found null walker so create the path in (relative to) pwd 
            if [ "x$WALKER" = "x" ]; then
                INVALID="ok"
                WALKER="DONE"
            fi
        fi
    done
    
    #   dbl check for accidentally writing over another installation
    #   skip this test if the path is not writeable
    if [ "$INVALID" = "ok" -a -d "$INSTALLBASE" ]; then
        while [ "$ANS" != "y" -a "$ANS" != "Y" -a "$ANS" != "n" -a "$ANS" != "N" ]
        do
            if [ $DONOTPROMPT = 0 ]
            then
              echo "WARNING:  $INSTALLBASE already exists.  Some files may be overwritten."
              echo "WARNING:  Do you wish to continue?  [y/n]"
              read ANS
            else
              ANS="Y"
            fi
        done
        if [ "$ANS" = "N" -o "$ANS" = "n" ]; then
            INVALID="INVALID"
            INSTALLBASE="?"
        fi
    fi 
    if [ $INVALID = "INVALID" -a $DONOTPROMPT = 1 ]
    then
       exit 1
    fi  
done

echo "Xyce will be installed in $INSTALLBASE"

echo "Decompressing archive and copying files..."
mkdir -p "$INSTALLBASE"
if [ "$?" -ne "0" ]; then
    echo "ERROR:  There was a problem creating $INSTALLBASE.  Exiting." 
    exit 1
fi

cd "$INSTALLBASE"

INSTALLBASE=`pwd`

cp -f "$ORIG_PWD/xtl.tar.gz" "$INSTALLBASE"

if [ "$?" -ne "0" ]; then
    echo "ERROR:  There was a problem copying files to the installation directory.  Exiting." 
    exit 1
fi

gzip -d "$INSTALLBASE/xtl.tar.gz"
if [ "$?" -ne "0" ]; then 
    echo "ERROR:  There was a problem decompressing the archive.  Exiting."
    exit 1
fi

tar xf "$INSTALLBASE/xtl.tar" > /dev/null 2>&1
if [ "$?" -ne "0" ]; then
    echo "ERROR:  There was a problem extracting files from the archive.  Exiting." 
    exit 1
fi


#   clean up (note that if user tries to install in $ORIG_PWD the tarball
#   gets wiped.  user can recover by unpacking installer again
rm -f "$INSTALLBASE/xtl.tar"

echo "installing Xyce here:   $INSTALLBASE" >> "$CONFIG_NOTES"


cd "$ORIG_PWD"
echo "$CONFIG_NOTES contains configuration notes."
echo
echo "Xyce installation completed successfully."
echo 

exit 0
