import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('inputFilename', 'list', #HTauTauAnalysis_1_1_Sl2.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input file name"
)


options.parseArguments()


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.TFileService=cms.Service("TFileService",fileName=cms.string("OUT_"+options.inputFilename+".root"))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")

process.demo = cms.EDAnalyzer('SignalStudy',
                                 InputFile = cms.string(options.inputFilename+".txt"),
)


process.p = cms.Path(process.demo)
