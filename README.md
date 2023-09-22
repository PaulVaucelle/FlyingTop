# FlyingTop
//-----------------------------Paul Vaucelle----------------------------//

//   -  Ph.D student for CMS, @IPHC, 10/2022->10/2025    -   //

//   ----------------           DisplacedTop analysis code      -----------------           //

// contact me @: paul.vaucelle@iphc.cnrs.fr -------------- //

// --------------- or paul.vaucelle@cern.ch      ---------------------      //

// ---------------- Supervisor : Daniel Bloch  --------------------   //

// ----------- contact @: <daniel.bloch@iphc.cnrs.fr>; -------- //

//----------------------------------------------------------------------------//

This LLP analysis code is available and can run on MiniAOD data tier.

For people @IPHC that want access to the code:

git clone https://github.com/Threshic/FlyingTop.git will only give these files and not the CMSSW env

Fro Run 2 Analyses, the release 10_6_30 is recommended. The analysis code is available for 10_6_20 (master branch)
and 10_6_30 (10_6_30 branch)

The following commands are for 10_6_20 but one can replace it by 10_6_30 as it is only installing the release.

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

mkdir FlyingTop
cd FlyingTop
mkedanlzr FlyingTop

git clone https://github.com/PaulVaucelle/FlyingTop.git


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

###You have to create an outputroot directory in the test directory if not already done

root

.L TreeReader.C+

### If the message "creating library " appears, you can try :

.L runTreeReader.C

### Else (it somehow happens that the command above crashes) copy paste the macro in the root terminal :)






