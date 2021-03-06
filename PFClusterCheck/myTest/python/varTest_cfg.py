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
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/14AEC6C4-D0C3-E411-BABE-0025905A6088.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/2AFA7B9B-C9C3-E411-9824-002590596490.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/2E5ADE85-C8C3-E411-B136-0025905A60F4.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/50DBC05C-C6C3-E411-98BC-0025905A6118.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/6EE51C37-D6C3-E411-AD4F-0026189438DF.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/96756E99-C9C3-E411-B059-0026189438FA.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/C8CDEDD2-D6C3-E411-9AC7-0026189438DF.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/DA4C090A-C6C3-E411-A285-0025905B858C.root',
        
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/1AEA54E8-C8C3-E411-BE74-0025905A6132.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/363AB411-D9C3-E411-99AF-0026189438E4.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/743E6681-D9C3-E411-AD6E-0025905B8596.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/8A322975-C8C3-E411-9D92-0025905A60CE.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/90E03F77-CEC3-E411-A6F7-0025905A60B8.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/94A9881A-CFC3-E411-9F8B-00261894391C.root',
        #'/store/relval/CMSSW_7_3_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/A8DB51BB-CEC3-E411-9AC0-0025905B858E.root',
        
        #'/store/relval/CMSSW_7_3_3/RelValSingleElectronPt10_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/56DAA656-C3C3-E411-9DD5-0025905B8596.root',

        #'/store/relval/CMSSW_7_3_3/RelValSingleElectronPt35_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/9EBEC239-CAC3-E411-B563-0025905A60EE.root',
        #'/store/relval/CMSSW_7_3_3/RelValSingleElectronPt35_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/A8C38233-CAC3-E411-90E5-003048FFD796.root',


        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/0272D672-CAC3-E411-B33B-002618943916.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/263676D0-D5C3-E411-A91D-0025905938B4.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/28F1E190-C9C3-E411-BC99-0025905B858A.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/343E05DA-C7C3-E411-9209-002590596490.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/467B0AC8-D5C3-E411-BC20-0025905A610A.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/AE71AD9B-C9C3-E411-A342-0025905964A6.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/B4D30004-C9C3-E411-8683-0025905A48EC.root',
        #'/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU50ns_MCRUN2_73_V10-v1/00000/DAF3C109-C6C3-E411-B485-0025905A6118.root',

#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/20B59843-D1C3-E411-9F09-0025905A60F4.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/3A8C2170-D0C3-E411-84B5-0025905A6104.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/3EEDFCAC-CBC3-E411-90E3-002618943957.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/42828BE3-E0C3-E411-ACFC-002618943925.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/46546A43-D1C3-E411-9931-0025905A60AA.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/60C0B0ED-E1C3-E411-A007-0025905B85E8.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/900A2BA6-CBC3-E411-9A29-0025905964C2.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/D6229421-CFC3-E411-90EB-0025905A612E.root',
#        '/store/relval/CMSSW_7_3_3/RelValH130GGgluonfusion_13/GEN-SIM-RECO/PU25ns_MCRUN2_73_V11-v1/00000/E8465DD1-CFC3-E411-A223-002618943826.root',
        #'/store/relval/CMSSW_7_3_3/RelValSingleGammaPt10_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/723E9440-C6C3-E411-A9A4-0025905A6110.root',
        #'/store/relval/CMSSW_7_3_3/RelValSingleGammaPt10_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/C48B5F3A-C6C3-E411-9C34-002590596490.root',

        #'/store/relval/CMSSW_7_3_1_patch1/RelValSingleGammaPt10_UP15/GEN-SIM-RECO//MCRUN2_73_V9_GenSim_7113-v1/00000/365463A4-DBB0-E411-BC39-0025905A609E.root',
        #'/store/relval/CMSSW_7_3_1_patch1/RelValSingleGammaPt10_UP15/GEN-SIM-RECO//MCRUN2_73_V9_GenSim_7113-v1/00000/7CF8ABA8-DBB0-E411-9672-0025905A606A.root',

        #'/store/relval/CMSSW_7_3_3/RelValSingleGammaPt35_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/52C10020-CCC3-E411-8BB1-0025905A60B8.root',
        #'/store/relval/CMSSW_7_3_3/RelValSingleGammaPt35_UP15/GEN-SIM-RECO/MCRUN2_73_V11-v1/00000/8E76F11F-CCC3-E411-AC23-0025905A60B8.root',


'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/0206BDC3-44DD-E411-B6CD-0025905B858C.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/02A71498-3BDD-E411-BC70-0025905938A4.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/0683A845-3DDD-E411-B1B5-0025905B8582.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/16F932AA-3BDD-E411-88F4-0025905B855C.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/185760C4-3EDD-E411-B4D0-0025905B85B2.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/2A067A44-3EDD-E411-93DA-003048FFCB8C.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/509433CB-44DD-E411-9B21-0025905B8606.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/56A4A49A-39DD-E411-8BD7-0025905A6090.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/58883EEF-3CDD-E411-8A2D-00261894385D.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/A6EA13B7-3ADD-E411-9A38-0025905938AA.root',
'/store/relval/CMSSW_7_4_0/RelValQCD_FlatPt_15_3000HS_13/GEN-SIM-RECO/MCRUN2_74_V7_gensim_740pre7-v1/00000/C2062ECB-38DD-E411-8654-003048FF9AA6.root',
   )
)

process.demo = cms.EDAnalyzer('varTest',
                              OutputFileName = cms.string("varTestQCD_740.root"),
                              label = cms.InputTag("gedPhotons"),
                              recHitsEBLabel = cms.InputTag('ecalRecHit', 'EcalRecHitsEB'),
                              recHitsEELabel = cms.InputTag('ecalRecHit', 'EcalRecHitsEE'),
)

process.p = cms.Path(process.demo)
