#!/bin/bash
# F.Poitevin - SLAC - June 2019
#
if [ ! -f ~/.osfrc ]; then 
  echo "Warning! ~/.osfrc was not found. This might prevent osf from working..."
else
  source ~/.osfrc; fi
echo "----------"; echo "OSF CLIENT"; echo "----------"
#
if [ $# -lt 1 ]; then
  echo "Usage error: $0 <local directory> [optional: <update list (yes/no)> <osf list file>]"; exit; fi
localdir=${1%/}
if [ ! -d $localdir ]; then
  echo "Usage error: $localdir does not exists."; exit; fi
updatels=${2:-'no'}
osfls=${3:-'osf_list.txt'}
if [ ! -f $osfls ]; then
  osf ls > $osfls
else
  if [ "$updatels" == 'yes' ]; then
    rm -f $osfls
    osf ls > $osfls
  fi
fi
echo "head $osfls : "
head $osfls
#
### Scenario 1: directory does not exist in OSF
if ! grep -q "$localdir" $osfls; then
  echo "$localdir is being created in OSF..."
  osf upload -r $localdir /
### Scenario 2: directory exists. Let's check it's complete
else
  echo "uploading to $localdir in OSF..."
  for file in $localdir/*; do
    if ! grep -q "$file" $osfls; then
      echo ">>> $file ..."
      osf upload "$file" "/$file"; fi
  done
fi
#
echo "----------"; echo "OSF CLIENT"; echo "----------"
