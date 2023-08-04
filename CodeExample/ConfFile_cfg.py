import FWCore.ParameterSet.Config as cms


from Configuration.AlCa.GlobalTag import GlobalTag # ????????
#from Configuration.StandardSequences.CondDBESSource_cff import GlobalTag # ??????

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process = cms.Process("Electrons")


dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)


process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames =
#cms.untracked.vstring('/store/mc/RunIISummer19UL18MiniAODv2/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/TrkExtra_106X_upgrade2018_realistic_v16_L1v1-v1/260000/0B9BA3E8-CA41-3742-ABB7-6AB02DFF5C03.root') #MC JPsi
cms.untracked.vstring('/store/data/Run2022F/ParkingSingleMuon0/MINIAOD/PromptReco-v1/000/360/390/00000/2c1a8864-bc25-4b39-9601-b2d6250652b4.root') #data BParking
#cms.untracked.vstring('/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/360/390/00000/05607e2d-f8e5-41a7-8392-bb5d44c16a6d.root') #data BParking
)

# define which IDs we want to produce
my_id_modules = [
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_iso_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_noIso_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff'
]

for idmod in my_id_modules:
   setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)

# The Ids produced are then accessible as edm::ValueMap<bool> objects, which can be accessed using the following input tags.

valueMapTags = {
  "eleVeto": cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-veto"),
  "eleLoose": cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-loose"),
  "eleMedium": cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-medium"),
  "eleTight": cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V2-tight"),
  "eleIsoWp80": cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp80"),
  "eleIsoWp90": cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp90"),
  "eleNoIsoWp80": cms.InputTag("mvaEleID-RunIIIWinter22-noIso-V1-wp80"),
  "eleNoIsoWp90": cms.InputTag("mvaEleID-RunIIIWinter22-noIso-V1-wp90")
}

process.electrons = cms.EDAnalyzer('ElectronAnalyzer',elecSrc = cms.untracked.InputTag("slimmedElectrons"),rhoSrc = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                                   pileupSrc = cms.untracked.InputTag("slimmedAddPileupInfo"),
                                   PV_Src = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),SV_Src = cms.untracked.InputTag("slimmedSecondaryVertices"),
                                   eleVeto= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-veto"),
                                   eleLoose= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-loose"),
                                   eleMedium= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V1-medium"),
                                   eleTight= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Winter22-122X-V2-tight"),
                                   eleIsoWp80= cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp80"),
                                   eleIsoWp90= cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-iso-V1-wp90"),
                                   eleNoIsoWp80= cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-noIso-V1-wp80"),
                                   eleNoIsoWp90= cms.InputTag("egmGsfElectronIDs:mvaEleID-RunIIIWinter22-noIso-V1-wp90")) # untracked.
                                   


process.TFileService = cms.Service("TFileService", fileName=cms.string("TnP_tree.root"))

##process.p = cms.Path(process.electrons*process.egammaPostRecoSeq)
process.p = cms.Path(process.electrons)


