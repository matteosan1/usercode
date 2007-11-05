#! /bin/bash

# 1 nome directory
# 2 numero eventi
# 3 process number

#mkdir /afs/cern.ch/user/s/sani/scratch1/CMSSW_1_5_2/src/NtupleMakers/ElectronAnalyzer/test/$1
cd /afs/cern.ch/user/s/sani/scratch1/CMSSW_1_6_0/src/NtupleMakers/ElectronAnalyzer/test/$1

export skip=$(( $3 * $2 ))

cat ../zeeC.cfg | sed -e "s/NEVENTS/${2}/g" | sed -e "s/SKIPEVENTS/${skip}/g" | sed -e "s/NOMEFILE/${1}_${3}/g" > ./${1}_${3}.cfg
eval `scramv1 runtime -sh`
cmsRun ${1}_${3}.cfg > output_${3}
