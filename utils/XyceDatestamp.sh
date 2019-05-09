#!/bin/sh

SRCDIR=$1
DATESTR=`date +%Y%m%d%H%M`
cd $SRCDIR
if [ -e .git ]
then
    GITSHA=`git log --pretty=format:%h -1`
    DATESTR="${DATESTR}-g-${GITSHA}"
fi
echo $DATESTR
