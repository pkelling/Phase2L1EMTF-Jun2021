import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process("Whatever",Phase2C9)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')  # load from DB
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
print "Using GlobalTag: %s" % process.GlobalTag.globaltag.value()


process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(321988),
    lastValue = cms.uint64(321988),
    interval = cms.uint64(1)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))


process.analyzer1 = cms.EDAnalyzer("MakePh2CoordLUT",
    # Verbosity level
    verbosity = cms.untracked.int32(1),
)

process.path1 = cms.Path(process.analyzer1)
