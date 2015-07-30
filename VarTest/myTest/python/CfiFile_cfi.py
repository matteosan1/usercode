import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('myTest',
                      OutputFileName = cms.string("eleOutput.root"),
                      label = cms.InputTag("slimmedElectrons")
)
