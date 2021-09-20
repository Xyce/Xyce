#!/bin/sh

SRCDIR=$1
DATESTR=`date +%Y%m%d%H%M`
cd $SRCDIR
if [ -e .git ]
then
    GITSHA=`git describe --dirty`
    DATESTR="${DATESTR}-(${GITSHA})"
fi
echo $DATESTR
