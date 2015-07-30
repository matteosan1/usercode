import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #'/store/relval/CMSSW_7_4_0/RelValSingleGammaPt10_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/384D6DF1-7ADD-E411-A663-0025905A60CA.root',

        '/store/relval/CMSSW_7_4_0/RelValSingleElectronPt1000_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/84880A63-30DD-E411-907E-0025905A606A.root',
        '/store/relval/CMSSW_7_4_0/RelValSingleElectronPt1000_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/F87DA85B-30DD-E411-B79F-0026189438F3.root',
        '/store/relval/CMSSW_7_4_0/RelValSingleElectronPt1000_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/FA885758-31DD-E411-A720-002618943896.root',
        

#'/store/relval/CMSSW_7_4_0/RelValSingleGammaPt35_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/3EECA996-26DD-E411-844F-0025905B857C.root',
#        '/store/relval/CMSSW_7_4_0/RelValSingleGammaPt35_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/5C986372-26DD-E411-995A-0025905A60B2.root',



        #'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt35_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/6E38A56F-26DD-E411-A617-0025905A60B6.root',
        #'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt35_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/94BF686D-26DD-E411-8D8B-0025905B8582.root',
        #'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt35_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7_gensim_740pre7-v1/00000/A62EDAFC-49DD-E411-B361-0025905B858E.root',


        #'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt10_UP15/GEN-SIM-DIGI-RAW-HLTDEBUG//MCRUN2_74_V7_gensim_740pre7-v1/00000/E415A7BA-26DD-E411-A70F-0025905B860C.root',

   )
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V7', '') #MCRUN2_73_V11', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi")
process.particleFlowClusterECAL.energyCorrector.applyMVACorrections = False

process.endjob_step = cms.EndPath(process.endOfProcess)

process.demo = cms.EDAnalyzer('myTest',
                              OutputFileName = cms.string("varTestEle1000_740_apr22_oldCorr.root"),
                              #label = cms.InputTag("gedPhotons"),
                              #isElectron = cms.bool(False),
                              label = cms.InputTag("gedGsfElectrons"),
                              isElectron = cms.bool(True),
)


process.p = cms.Path(process.demo)

from SLHCUpgradeSimulations.Configuration.postLS1Customs import *
process = customisePostLS1(process)
