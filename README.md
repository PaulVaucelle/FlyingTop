FlyingTop

//-----------------------------Paul Vaucelle----------------------------//

// - Ph.D student for CMS, @IPHC, 10/2022->10/2025 - //

// ---------------- DisplacedTop analysis code ----------------- //

// contact me @: paul.vaucelle@iphc.cnrs.fr -------------- //

// --------------- or paul.vaucelle@cern.ch --------------------- //

// ---------------- Supervisor : Daniel Bloch -------------------- //

// ----------- contact @: daniel.bloch@iphc.cnrs.fr; -------- //

//----------------------------------------------------------------------------//

The analysis code is available and can run on MiniAOD data tier. For people @IPHC that want access to the code:

git clone https://github.com/Threshic/FlyingTop.git will only give these files and not the CMSSW env

For Run 3 Analyses, we are currently ousing the release 14_0_8. The analysis code is available for Run 2 : 10_6_30 (10_6_30 is the default branch) and Run 3 : 14_0_8
Release 14_0_8_FLY :
create repository :

FIrst time :

// Dans un envrionnement CMS (10_6_30):
export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700

eval `scramv1 runtime -sh`

export PYTHONPATH=$PYTHONPATH:/opt/glite/lib64/python:/opt/fpconst/lib/python2.4/site-packages
export PYTHONPATH=${PYTHONPATH}:${GLITE_LOCATION}/lib64
cmsenv

export SCRAM_ARCH=el8_amd64_gcc12
 // scram list -a # pour voir les releases si cela t'intéresse
cmssw-el8 # ouvre la singularité
// revenir à l'endroit où tu veux mettre la nouvelle release
scramv1 p CMSSW CMSSW_14_0_8
mv CMSSW_14_0_8 CMSSW_14_0_8_FLY
cd CMSSW_14_0_8_FLY
scramv1 b ProjectRename
cd src
eval  `scramv1 r -sh`

git cms-init

git cms-addpkg RecoEgamma/EgammaTools

git clone https://github.com/cms-egamma/EgammaPostRecoTools.git

mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.

git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/

git cms-addpkg EgammaAnalysis/ElectronTools

git cms-addpkg RecoJets/JetProducers

scramv1 b clean
scramv1 b -j4


// Par la suite :


cd CMSSW_14_0_8/src
export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el8_amd64_gcc12
cmssw-el8
eval `scramv1 runtime -sh`
export PYTHONPATH=$PYTHONPATH:/opt/glite/lib64/python:/opt/fpconst/lib/python2.4/site-packages
export PYTHONPATH=${PYTHONPATH}:${GLITE_LOCATION}/lib64
// après ce point, on peut rentrer les commandes habituellles:
scramv1 b clean
scramv1 b -j4

Electron p4 and Muon rochester corrections are directly applied in the code through different methods as well as for PileUp
compilation :

scramv1 b clean

scramv1 b -j6


source /cvmfs/cms.cern.ch/crab3/crab.sh voms-proxy-init -rfc -voms cms -valid 192:00

cd FlyingTop/FlyingTop/test

cmsRun flyingtop_MC.py or cmsRun flyingtop_data.py (specific to the branch 10_6_30)

or (still working but not used, see next point "To run jobs")

crab submit -c crab_config_mc_2018.py (use grid : don't forget the certificate)
To run jobs :

./multicrab --crabCmd submit (edit to choose the samples)

./multicrab --crabCmd status --workArea ./<work_directory> (work directory is the date of the launching of the jobs DD_MM_YYYY)

./multicrab --crabCmd getoutput --workArea ./<work_directory> --crabCmdOpts --checksum=no (to retrieve the root files)
To merge the outputfiles, use the haddWithWeights.py macro:

Please check the paths so that they correspond to your envrionment

python haddWithWeights.py <work_directory>/

