import FWCore.ParameterSet.Config as cms

shallowGainCalibration = cms.EDProducer("ShallowGainCalibration",
                                      Tracks=cms.InputTag("generalTracks",""),
                                      Clusters=cms.InputTag("siStripClusters"),
                                      Prefix=cms.string("GainCalibration"),
                                      Suffix=cms.string(""))
