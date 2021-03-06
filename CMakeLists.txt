# Copyright CERN and copyright holders of ALICE O2. This software is distributed
# under the terms of the GNU General Public License v3 (GPL Version 3), copied
# verbatim in the file "COPYING".
#
# See http://alice-o2.web.cern.ch/license for full licensing information.
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.

# FIXME Is this one supposed to be a header only library (in which case the .h
# to be installed should be in include/TestWorkflows) or not a library at all ?
# o2_add_library(TestWorkflows SOURCES src/dummy.cxx
#
# set(HEADERS src/o2_sim_its_ALP3.h src/o2_sim_tpc.h )
#

o2_add_executable(track-collection-iteration
                  SOURCES src/trackCollectionIteration.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(track-iteration
                  SOURCES src/trackIteration.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(full-track-iteration
                  SOURCES src/fullTrackIteration.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(collision-tracks-iteration
                  SOURCES src/collisionTracksIteration.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(new-collections
                  SOURCES src/newCollections.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(associated-derived
                  SOURCES src/associatedExample.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(dynamic-columns
                  SOURCES src/dynamicColumns.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(histograms
                  SOURCES src/histograms.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(filters
                  SOURCES src/filters.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                  COMPONENT_NAME AnalysisTutorial)

o2_add_executable(efficiency-analysis
                  SOURCES src/efficiencyAnalysis.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                                        O2::SimulationDataFormat
                                        O2::ReconstructionDataFormats
                                        O2::GlobalTracking
                                        O2::DataFormatsTOF
                                        O2::TOFSimulation
                  COMPONENT_NAME AnalysisTutorial)
