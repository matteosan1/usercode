import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("HLTTest.MyCandidates.turnon_cfi")
process.turnOn.RootFileName = cms.string("hltD.root")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/data/Run2012D/DoubleElectron/RAW-RECO/ZElectron-PromptSkim-v1/000/207/886/00000/300F113D-E238-E211-BBD9-002618FDA210.root',
    )
                            )

process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_53_V21::All"
#process.GlobalTag.globaltag = "GR_P_V42::All"

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.p1 = cms.Path(process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)
process.p2 = cms.Path(process.turnOn)
