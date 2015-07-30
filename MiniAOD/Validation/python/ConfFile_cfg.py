import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8_unsch-v1/00000/1A9854B1-C806-E511-BB24-003048FFD752.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8_unsch-v1/00000/1E89ECEB-CE06-E511-B047-0025905A60E0.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8_unsch-v1/00000/2A0D05D2-3507-E511-A2B9-0025905964C4.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8_unsch-v1/00000/4A31CFAA-3507-E511-A20B-0025905A48E4.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9_unsch-v1/00000/16E9F3F9-3307-E511-8D85-0025905A60BE.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9_unsch-v1/00000/90BE2710-DC06-E511-B752-0025905B85D8.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9_unsch-v1/00000/9472E50A-3407-E511-A008-0025905964C4.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9_unsch-v1/00000/D63DFDFD-DA06-E511-B401-00261894393A.root',
        
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/MCRUN2_74_V9_unsch-v1/00000/0C156DCB-D106-E511-A508-0025905A611C.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/MCRUN2_74_V9_unsch-v1/00000/348884C8-D106-E511-A058-0025905A60B0.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/GEN-SIM-RECO/MCRUN2_74_V9_unsch-v1/00000/D454C277-CD06-E511-9A52-0025905A6090.root',
        
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_unsch-v1/00000/A8235305-3407-E511-8ED9-0025905A60B4.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_unsch-v1/00000/E826F00B-3407-E511-9B7A-0025905A612C.root',

        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9_unsch-v1/00000/B423E53C-D206-E511-B6AB-0025905B85E8.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9_unsch-v1/00000/C6B18245-D206-E511-ADEC-0025905A60EE.root',

        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/MINIAODSIM/PU50ns_MCRUN2_74_V8_unsch-v1/00000/2EB0AC8D-3507-E511-AF76-003048FFCB9E.root',
        #'/store/relval/CMSSW_7_4_3_patch1/RelValZEE_13/MINIAODSIM/PU50ns_MCRUN2_74_V8_unsch-v1/00000/B4C1527D-3507-E511-A0D8-0025905964C4.root',
        

        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9-v9/00000/8C84085F-1306-E511-9438-0025905A60BC.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/MCRUN2_74_V9-v9/00000/E0371A55-1306-E511-A9E5-0025905964C0.root',
        
        
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PU25ns_MCRUN2_74_V9-v11/00000/28E326F3-5A07-E511-9AE4-003048FF9AA6.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PU25ns_MCRUN2_74_V9-v11/00000/86D08403-5B07-E511-9673-0025905A60CA.root',
        
        
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PU50ns_MCRUN2_74_V8-v11/00000/B81FC72C-3C07-E511-B2BA-0025905A607A.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PU50ns_MCRUN2_74_V8-v11/00000/D25E113E-3C07-E511-832F-0025905B85EE.root',
      
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/MCRUN2_74_V9-v9/00000/3068E195-0D06-E511-835B-0025905964A2.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/MCRUN2_74_V9-v9/00000/C84EEA95-0D06-E511-A3FE-0025905A610C.root',
        
        
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9-v11/00000/100907D1-E406-E511-BB65-00261894398A.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9-v11/00000/5C7CEA03-2307-E511-8FDE-0025905B8606.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9-v11/00000/985DF8AC-DB06-E511-B50E-0025905A610C.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9-v11/00000/A83C5DBC-5407-E511-8B92-0025905A6084.root',
        #'/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU25ns_MCRUN2_74_V9-v11/00000/A8810CB9-5407-E511-80AE-0025905A608A.root',
        
        
        '/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8-v11/00000/309E4BC9-3807-E511-9B2B-0025905A48D0.root',
        '/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8-v11/00000/685B87C6-D306-E511-9F00-0025905A497A.root',
        '/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8-v11/00000/7035570D-2107-E511-ACC0-003048FFCC0A.root',
        '/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8-v11/00000/981C953F-CF06-E511-84BD-003048FFD79C.root',
        '/store/relval/CMSSW_7_4_3/RelValZEE_13/GEN-SIM-RECO/PU50ns_MCRUN2_74_V8-v11/00000/9A612ED7-2407-E511-B089-0025905B85AA.root',



        )
)

process.demo = cms.EDAnalyzer('Validation',
                              #isMini = cms.bool(False),
                              outputFileName = cms.string("test.root"),
                              electrons = cms.InputTag("gedGsfElectrons")
                              #isMini = cms.bool(True),
                              #outputFileName = cms.string("miniaodsim_PU50ns_MCRUN2_74_V8-v11_sche.root")
)


process.p = cms.Path(process.demo)
