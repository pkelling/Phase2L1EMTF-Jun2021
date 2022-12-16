import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C12_cff import Phase2C12
process = cms.Process("GenEMTFLuts",Phase2C12)

# Use standard seq to load geom (?)
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')  # load from DB

# load geometries directly
process.load('Configuration.Geometry.GeometryExtended2026D70_cff')
process.load('Configuration.Geometry.GeometryExtended2026D70Reco_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
print("Using GlobalTag: ", process.GlobalTag.globaltag.value())

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))


process.analyzer1 = cms.EDAnalyzer("MakePh2CoordLUT",
    # Verbosity level
    verbosity = cms.untracked.int32(1),
)

process.source = cms.Source("EmptySource")

process.path1 = cms.Path(process.analyzer1)
