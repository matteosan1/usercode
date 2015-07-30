import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardConfiguration.Reconstruction_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.Digi_cff')
#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
#process.load('Configuration.StandardSequences.DigiToRaw_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V7', '') #MCRUN2_73_V11', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/0206BDC3-44DD-E411-B6CD-0025905B858C.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/02A71498-3BDD-E411-BC70-0025905938A4.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/0683A845-3DDD-E411-B1B5-0025905B8582.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/16F932AA-3BDD-E411-88F4-0025905B855C.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/185760C4-3EDD-E411-B4D0-0025905B85B2.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/2A067A44-3EDD-E411-93DA-003048FFCB8C.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/509433CB-44DD-E411-9B21-0025905B8606.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/56A4A49A-39DD-E411-8BD7-0025905A6090.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/58883EEF-3CDD-E411-8A2D-00261894385D.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/A6EA13B7-3ADD-E411-9A38-0025905938AA.root',
#'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/C2062ECB-38DD-E411-8654-003048FF9AA6.root',


#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/0277FEE7-45DD-E411-88BB-0025905B85AA.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/3019B174-45DD-E411-A1E5-003048FF9AC6.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/A4000EE9-45DD-E411-9D61-0025905938AA.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/AE64A931-4EDD-E411-90C1-0025905A6064.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/D429293B-48DD-E411-9863-0025905B85AA.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/D43CFF28-4EDD-E411-9022-0025905B855C.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/E28AEC3E-44DD-E411-9D5E-0025905B85AA.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/E4BE6A3F-44DD-E411-BCAC-0025905A6076.root',
#'/store/relval/CMSSW_7_4_0/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V7_gensim_740pre7-v1/00000/EE01239F-48DD-E411-9E02-002618943962.root',


#'/store/relval/CMSSW_7_4_0/RelValSingleGammaPt35_UP15/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/8020DDA3-5BDD-E411-A6DD-0025905A605E.root',
#'/store/relval/CMSSW_7_4_0/RelValSingleGammaPt35_UP15/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/9A0469A1-5BDD-E411-BD67-0026189438F9.root',

#'/store/relval/CMSSW_7_4_0/RelValSingleGammaPt10_UP15/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/0A4AECE3-D5DD-E411-B8C2-0025905B8572.root'

#'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt10_UP15/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/346CFC80-91DD-E411-A407-002590596484.root'

'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt35_UP15/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/46EA0189-76DD-E411-BDDD-0025905A6076.root',
'/store/relval/CMSSW_7_4_0/RelValSingleElectronPt35_UP15/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/621C5A8D-76DD-E411-BBA3-0025905964A6.root',


   )
)
process.load('RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff')
process.reclustering_step = cms.Path(process.pfClusteringPS * process.pfClusteringECAL)

process.demo = cms.EDAnalyzer('varTest',
                              OutputFileName = cms.string("varTestEle35_740_apr22.root"),
                              label = cms.InputTag("gedPhotons"),
                              eetops = cms.InputTag("particleFlowClusterECAL"),#,"PFClusterAssociationES"),
                              recHitsEBLabel = cms.InputTag('ecalRecHit', 'EcalRecHitsEB'),
                              recHitsEELabel = cms.InputTag('ecalRecHit', 'EcalRecHitsEE'),
)

process.p = cms.Path(process.demo)
