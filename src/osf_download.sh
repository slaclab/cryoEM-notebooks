#!/bin/bash
# F.Poitevin - SLAC - June 2019
#
# This script downloads the content of an OSF directory into a local directory
#
if [ ! -f ~/.osfrc ]; then
  echo "Warning! ~/.osfrc was not found. This might prevent osf from working..."
else
  source ~/.osfrc; fi
echo "----------"; echo "OSF CLIENT"; echo "----------"
#
if [ $# -lt 2 ]; then
  echo "Usage error: $0 <remote directory (keyword-based search)> <local directory (exact relative path)> [optional: <update list (yes/no)> <osf list file>]"; exit; fi
remotedir=${1%/}
localdir=$2
if [ ! -d $localdir ]; then
  echo "Please create local directory before running this..."
  echo "$ mkdir $localdir"
  exit
fi
updatels=${3:-'no'}
osfls=${4:-'osf_list.txt'}
if [ ! -f $osfls ]; then
  osf ls > $osfls
else
  if [ "$updatels" == 'yes' ]; then
    rm -f $osfls
    osf ls > $osfls ;fi; fi
echo "head $osfls : "
head $osfls
#
### Scenario 1: directory does not exist in OSF
if ! grep -q "$remotedir" $osfls; then
  echo "$remotedir was not found in OSF..."
### Scenario 2: directory exists. Let's check it's complete
else
  echo "downloading content of $remotedir from OSF into current directory..."
  grep "$remotedir" $osfls > "tmp_${remotedir}.txt"
  while read -r line
  do
    if [ ! -f ${localdir}/${line##*/} ]; then
      echo ">>> $line ..."
      osf fetch $line ${localdir}/${line##*/}
    fi
  done < tmp_${remotedir}.txt
fi
#
echo "----------"; echo "OSF CLIENT"; echo "----------"
