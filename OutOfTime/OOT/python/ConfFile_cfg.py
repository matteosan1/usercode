import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-pythia6/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/044C97C0-B675-E311-B49C-0025901D4D76.root',
        )
                            )

process.demo = cms.EDAnalyzer('OOT'
)


process.p = cms.Path(process.demo)
