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
For people @IPHC that want access to the code:

git clone https://github.com/Threshic/FlyingTop.git will only give these files and not the CMSSW env

Fro Run 2 Analyses, the release 10_6_30 is recommended. The analysis code is available for 10_6_20 
and 10_6_30 (10_6_30  is the default one)


# Release 10_6_30_FLY : 

###
### create repository :
###

cd /opt/sbg/cms/insert uiX_dataY/insert name

export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch

source /cvmfs/cms.cern.ch/cmsset_default.sh

source /cvmfs/cms.cern.ch/crab3/crab.sh 

export SCRAM_ARCH=slc7_amd64_gcc700

scramv1 p CMSSW CMSSW_10_6_30

mv CMSSW_10_6_30 CMSSW_10_6_30_FLY

cd CMSSW_10_6_30_FLY

scramv1 b ProjectRename

cd src

eval  `scramv1 r -sh`

git cms-init

mkdir FlyingTop

cd FlyingTop

mkedanlzr FlyingTop

git clone https://github.com/PaulVaucelle/FlyingTop.git

git checkout CMSSW_10_6_30

----------------------------

# packages to add :
# Do these steps  in CMSSW_10_6_30/src

cd CMSSW_10_6_30/src

git cms-addpkg RecoEgamma/EgammaTools

git clone https://github.com/cms-egamma/EgammaPostRecoTools.git

 mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
 
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/

git cms-addpkg EgammaAnalysis/ElectronTools

git cms-addpkg RecoJets/JetProducers

scram b -j 8

More details are available at : https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018 and https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL#Recommendations_for_2018_UL_data if needed
 

# these packages are used in the flyingtop_MC(data).py files 

# Electron p4 and Muon rochester corrections are directly applied in the code through different methods as well as for PileUp

###
### compilation :
###

scramv1 b clean

scramv1 b -j6 

-----------------------------------

###
### Execute :
###

cd /opt/sbg/cms/<insert uiX_dataY>/<insert name>/CMSSW_10_6_30_FLY/src

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

cmsRun flyingtop_MC.py or cmsRun flyingtop_data.py (specific to the branch 10_6_30)
  
 or (still working but not used, see next point "To run jobs")
  
 crab submit -c crab_config_mc_2018.py (use grid : don't forget the certificate)
 
# To run jobs :

./multicrab --crabCmd submit  (edit to choose the samples)

./multicrab --crabCmd status --workArea ./<work_directory>   (work directory is the date of the launching of the jobs DD_MM_YYYY)

./multicrab --crabCmd getoutput --workArea ./<work_directory> --crabCmdOpts --checksum=no  (to retrieve the root files)

# To merge the outputfiles, use the haddWithWeights.py macro:

Please check the paths so that they correspond to your envrionment

python haddWithWeights.py <work_directory>/

Please do not forget the "/" after <work_directory>

# To run TreeReader.C (current version gives efficiencies for the different selections and the reconstruction of vertices + a txt file to summarize it)

### You have to create an outputroot directory in the test directory if not already done

root

.L TreeReader.C+

### If the message "creating library " appears, you can try :

.L runTreeReader.C

### Else (it somehow happens that the command above crashes) copy paste the macro in the root terminal :)






