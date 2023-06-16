# FlyingTop
//-----------------------------Paul Vaucelle----------------------------//

//   -  Ph.D student for CMS, @IPHC, 10/2022->10/2025    -   //

//   ----------------           DisplacedTop analysis code      -----------------           //

// contact me @: paul.vaucelle@iphc.cnrs.fr -------------- //

// --------------- or paul.vaucelle@cern.ch      ---------------------      //

// ---------------- Supervisor : Daniel Bloch  --------------------   //

// ----------- contact @: <daniel.bloch@iphc.cnrs.fr>; -------- //

//----------------------------------------------------------------------------//

The analysis code is available and can run on MiniAOD data tier.

A Reco version is available @: https://github.com/Threshic/TrackingPerf

Look for the _RECO file or "Version-05_10_2022-AVF (Pour présentation TPS)" of TrackingPerf.cc (in plugins) on the history version of the repo.

For the _RECO file, there are two versions : 16/11/22 and 17/11/22. Both are quite the same, the first is more of a draft than the second one. You'd rather use the second one (Daniel).

For people @IPHC that want access to the code:

git clone https://github.com/Threshic/FlyingTop.git will only give these files and not the CMSSW env

# Release 10_6_20_FLY : 

###
### create repository :
###

cd /opt/sbg/cms/insert uiX_dataY/insert name

export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch

source /cvmfs/cms.cern.ch/cmsset_default.sh

source /cvmfs/cms.cern.ch/crab3/crab.sh 

export SCRAM_ARCH=slc7_amd64_gcc700

scramv1 p CMSSW CMSSW_10_6_20

mv CMSSW_10_6_20 CMSSW_10_6_20_FLY

cd CMSSW_10_6_20_FLY

scramv1 b ProjectRename

cd src

eval  `scramv1 r -sh`

git cms-init

 cp -r /opt/sbg/cms/ui2_data1/blochd/CMSSW_10_6_20_FLY_pourPaul/src/FlyingTop . (Daniel's version)
 
 or
 
cp -r /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop . (Paul's version)

just for information for the very first time:

mkdir FlyingTop
cd FlyingTop
mkedanlzr FlyingTop
edit => FlyingTop/plugins/BuildFile.xml & // and add :

use name="DataFormats/PatCandidates"/

###
### compilation :
###

scramv1 b clean

scramv1 b -j6 

-----------------------------------

###
### Execute :
###

cd /opt/sbg/cms/<insert uiX_dataY>/<insert name>/CMSSW_10_6_20_FLY/src

export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700

eval `scramv1 runtime -sh`

export PYTHONPATH=$PYTHONPATH:/opt/glite/lib64/python:/opt/fpconst/lib/python2.4/site-packages
export PYTHONPATH=${PYTHONPATH}:${GLITE_LOCATION}/lib64
cmsenv

# to specify crab version :
source /cvmfs/cms.cern.ch/crab3/crab.sh 
voms-proxy-init -rfc -voms cms -valid 192:00

cd FlyingTop/FlyingTop/test
cmsRun flyingtop.py (to run locally)
  
 or
  
 crab submit -c crab_config_mc_2018.py (use grid : don't forget the certificate)










