#!/bin/bash

#parse input arguments
YEAR=2017
while getopts "o:y:" opt; do
  case "$opt" in
	o) WHAT=$OPTARG
	  ;;
	y) YEAR=$OPTARG
	  ;;
  esac
done

#check an operation has been given
if [ -z "${WHAT}" ]; then
  echo "run.sh -o <LOCAL/SUBMIT> [ -y 2016/2017/2018 ] ";
  echo " LOCAL		- run locally";
  echo " SUBMIT 	- submit jobs";
  exit 1;
fi

case ${WHAT} in

  LOCAL)
	python $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/BH_postproc.py -m -i $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/TTWJetsToLNu_TuneCP5_13TeV_NanoAOD.root --year ${YEAR} -o $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test;
	
        #python $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/BH_postproc.py -e B -i $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/Run2017B_SingleMuon_NanoAODv2.root --year ${YEAR} -o $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test;

	#python $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/BH_postproc.py -e B -i $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/Run2017B_SingleMuon_NanoAODv2.root --year ${YEAR} -o $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test;

	#python $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/run.py -d -i $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test/Run2017B_SingleMuon_NanoAODv2.root --year ${YEAR} -o $CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/BHplus/test;
	;;

#  SUBMIT)
	
esac
