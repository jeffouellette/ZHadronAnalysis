#! /bin/bash

runFile=runMinbias18
echo "Universe = vanilla" > ${runFile}.job
echo "Notification = Complete" >> ${runFile}.job
echo "Executable = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Macros/ProcessFile" >> ${runFile}.job
echo "Priority = +1" >> ${runFile}.job
echo "GetEnv = True" >> ${runFile}.job
echo "InitialDir = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Macros" >> ${runFile}.job
echo "accounting_group = group_atlas.boulder" >> ${runFile}.job
echo "request_memory   = 8 GB" >> ${runFile}.job

echo "" >> ${runFile}.job

echo "sysFlags = false none" >> ${runFile}.job # options are "none", "doMixVar[A-H]", or "doPPMixVar"

echo "" >> ${runFile}.job
echo "" >> ${runFile}.job

for iRG in 'GroupA' 'GroupB' 'GroupC' 'GroupD' 'GroupE' 'GroupF' 'GroupG' 'GroupH' 'GroupI' 'GroupJ' 'GroupK'
do
  for iCent in {1..28}
  do
    echo "Arguments = minbias ${iRG}/${iRG}_iCent${iCent}.root 0 ${iRG}/${iRG}_iCent${iCent}_hists.root 0 true \$(sysFlags)" >> ${runFile}.job
    echo "Error = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/errors/${iRG}_iCent${iCent}.err" >> ${runFile}.job
    echo "Log = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/logs/${iRG}_iCent${iCent}.log" >> ${runFile}.job
    echo "Output = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/outputs/${iRG}_iCent${iCent}.out" >> ${runFile}.job
    echo "Queue 1" >> ${runFile}.job
    echo "" >> ${runFile}.job
  done
done



runFile=runMCMixing
echo "Universe = vanilla" > ${runFile}.job
echo "Notification = Complete" >> ${runFile}.job
echo "Executable = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Macros/ProcessFile" >> ${runFile}.job
echo "Priority = +1" >> ${runFile}.job
echo "GetEnv = True" >> ${runFile}.job
echo "InitialDir = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Macros" >> ${runFile}.job
echo "accounting_group = group_atlas.boulder" >> ${runFile}.job
echo "request_memory   = 8 GB" >> ${runFile}.job

echo "" >> ${runFile}.job

echo "sysFlags = false none" >> ${runFile}.job # options are "none", "doMixVar[A-H]", or "doPPMixVar"

echo "" >> ${runFile}.job
echo "" >> ${runFile}.job

echo "Arguments = mcminbias CombinedTrees/pp_all.root pp_all.root MixedHistograms/pp_all_hists.root Pythia/pp_all.root false \$(sysFlags)" >> ${runFile}.job
echo "Error = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/errors/pp_all.err" >> ${runFile}.job
echo "Log = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/logs/pp_all.log" >> ${runFile}.job
echo "Output = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/outputs/pp_all.out" >> ${runFile}.job
echo "Queue 1" >> ${runFile}.job
echo "" >> ${runFile}.job

for iRG in 'GroupD' 'GroupE' 'GroupF' 'GroupG' 'GroupH' 'GroupI' 'GroupJ' 'GroupK'
do
  for iCent in {1..28}
  do
    echo "Arguments = mcminbias CombinedTrees/${iRG}_iCent${iCent}.root ${iRG}/${iRG}_iCent${iCent}.root MixedHistograms/${iRG}_iCent${iCent}_hists.root DataOverlay/${iRG}/${iRG}_iCent${iCent}.root true \$(sysFlags)" >> ${runFile}.job
    echo "Error = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/errors/${iRG}_iCent${iCent}.err" >> ${runFile}.job
    echo "Log = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/logs/${iRG}_iCent${iCent}.log" >> ${runFile}.job
    echo "Output = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/outputs/${iRG}_iCent${iCent}.out" >> ${runFile}.job
    echo "Queue 1" >> ${runFile}.job
    echo "" >> ${runFile}.job
  done
done



runFile=runMC
echo "Universe = vanilla" > ${runFile}.job
echo "Notification = Complete" >> ${runFile}.job
echo "Executable = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Macros/ProcessFile" >> ${runFile}.job
echo "Priority = +1" >> ${runFile}.job
echo "GetEnv = True" >> ${runFile}.job
echo "InitialDir = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Macros" >> ${runFile}.job
echo "accounting_group = group_atlas.boulder" >> ${runFile}.job

echo "" >> ${runFile}.job

echo "sysFlags = false none" >> ${runFile}.job # options are "none", "doMixVar[A-H]", or "doPPMixVar"

echo "" >> ${runFile}.job
echo "" >> ${runFile}.job

echo "Arguments = mc CombinedTrees/pp_all.root 0 Histograms/pp_all_hists.root 0 false \$(sysFlags)" >> ${runFile}.job
echo "Error = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/errors/pp_all.err" >> ${runFile}.job
echo "Log = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/logs/pp_all.log" >> ${runFile}.job
echo "Output = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/outputs/pp_all.out" >> ${runFile}.job
echo "Queue 1" >> ${runFile}.job
echo "" >> ${runFile}.job

for iRG in 'GroupD' 'GroupE' 'GroupF' 'GroupG' 'GroupH' 'GroupI' 'GroupJ' 'GroupK'
do
  for iCent in {1..28}
  do
    echo "Arguments = mc CombinedTrees/${iRG}_iCent${iCent}.root 0 Histograms/${iRG}_iCent${iCent}_hists.root 0 true \$(sysFlags)" >> ${runFile}.job
    echo "Error = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/errors/${iRG}_iCent${iCent}.err" >> ${runFile}.job
    echo "Log = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/logs/${iRG}_iCent${iCent}.log" >> ${runFile}.job
    echo "Output = /atlasgpfs01/usatlas/workarea/jeff/atlas-hi/ZTrackAnalysis/Submit/outputs/${iRG}_iCent${iCent}.out" >> ${runFile}.job
    echo "Queue 1" >> ${runFile}.job
    echo "" >> ${runFile}.job
  done
done
