#!/bin/tcsh

eval `scramv1 runtime -csh`
rehash


echo " "
echo " Production chain (four jobs: digi+digi2raw, HLT, split, reco):"

echo " "
echo "/bin/rm RelVal_Digi_Raw.root RelVal_Pure_Raw.root RelVal_Digi_Digi2Raw.log"
      /bin/rm RelVal_Digi_Raw.root RelVal_Pure_Raw.root RelVal_Digi_Digi2Raw.log
echo "./testcfg RelVal_Digi_Digi2Raw.cfg   >& RelVal_Digi_Digi2Raw.log"
      ./testcfg RelVal_Digi_Digi2Raw.cfg   >& RelVal_Digi_Digi2Raw.log
#echo "cmsRun --strict RelVal_Digi_Digi2Raw.cfg   >& RelVal_Digi_Digi2Raw.log"
#      cmsRun --strict RelVal_Digi_Digi2Raw.cfg   >& RelVal_Digi_Digi2Raw.log

echo " "
echo "/bin/rm HLTFromDigiRaw.root             RelVal_HLTFromRaw.log"
      /bin/rm HLTFromDigiRaw.root             RelVal_HLTFromRaw.log
echo "./testcfg RelVal_HLTFromRaw.cfg      >& RelVal_HLTFromRaw.log"
      ./testcfg RelVal_HLTFromRaw.cfg      >& RelVal_HLTFromRaw.log
#echo "cmsRun --strict RelVal_HLTFromRaw.cfg      >& RelVal_HLTFromRaw.log"
#      cmsRun --strict RelVal_HLTFromRaw.cfg      >& RelVal_HLTFromRaw.log

echo " "
echo "/bin/rm CSA07*.root                     RelVal_PrimaryDatasets.log"
      /bin/rm CSA07*.root                     RelVal_PrimaryDatasets.log
echo "./testcfg RelVal_PrimaryDatasets.cfg >& RelVal_PrimaryDatasets.log"
      ./testcfg RelVal_PrimaryDatasets.cfg >& RelVal_PrimaryDatasets.log

echo " "
echo "/bin/rm RelVal_Reco.root                RelVal_Reco.log"
      /bin/rm RelVal_Reco.root                RelVal_Reco.log
echo "./testcfg RelVal_Reco.cfg            >& RelVal_Reco.log"
      ./testcfg RelVal_Reco.cfg            >& RelVal_Reco.log
#echo "cmsRun --strict RelVal_Reco.cfg            >& RelVal_Reco.log"
#      cmsRun --strict RelVal_Reco.cfg            >& RelVal_Reco.log


echo " "
echo " Quick test running HLTtable from pure Raw:"

echo " "
echo "/bin/rm HLTFromPureRaw.root       HLTtable.log"
      /bin/rm HLTFromPureRaw.root       HLTtable.log
echo "./testcfg HLTtable.cfg         >& HLTtable.log"
      ./testcfg HLTtable.cfg         >& HLTtable.log
#echo "cmsRun --strict HLTtable.cfg         >& HLTtable.log"
#      cmsRun --strict HLTtable.cfg         >& HLTtable.log

echo " "
echo "Finished!"
