import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

process = cms.Process("GenEMTFLuts",Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
print "Using GlobalTag: %s" % process.GlobalTag.globaltag.value()

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))


process.analyzer1 = cms.EDAnalyzer("MakePh2CoordLUT",
    # Verbosity level
    verbosity = cms.untracked.int32(1),
)

process.source = cms.Source("EmptySource")

process.path1 = cms.Path(process.analyzer1)
