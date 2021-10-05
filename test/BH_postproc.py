import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.BHplus.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.BHplus.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.BHplus.modules.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.BHplus.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.BHplus.modules.BHProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=False, action='store_true')
  parser.add_option('-e', dest='era',default="B", help="Run period, only work for data")
  parser.add_option('-i', '--in', dest='inputs', help='input directory with files', default=None, type='string')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  parser.add_option('-o', '--out', dest='output', help='output directory with files', default=None, type='string')

  (opt, args) = parser.parse_args()

  if opt.ismc:
    if opt.year == "2016" and opt.era == 'A':
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puAutoWeight_2016(),PrefCorr(),muonIDISOSF2016(),muonScaleRes2016a(),eleRECOSF2016(),eleIDSF2016(), BH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016" and opt.era == 'B':
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puAutoWeight_2016(),PrefCorr(),muonIDISOSF2016(),muonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(), BH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017":
      #p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puWeight_2017(),PrefCorr(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(), jmeCorrections_UL2017MC(),btagSF2017UL(), BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
      p = PostProcessor(".", [opt.inputs], modules=[countHistogramsModule(),puWeight_2017(),PrefCorr(),muonIDISOSF2017(),muonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(), jmeCorrections_UL2017MC(),btagSF2017UL(), BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt",maxEntries=1000)
    if opt.year == "2018":
      p = PostProcessor(".", inputFiles(), modules=[countHistogramsModule(),puAutoWeight_2018(),muonIDISOSF2018(),muonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(), BH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016" and (opt.era in ["B","C","D"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),BH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016" and (opt.era in ["E","F"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016a(),BH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016" and (opt.era in ["G","H"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2016b(),BH2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")

    if opt.year == "2017" and (opt.era in ["B"]):
      p = PostProcessor(".", [opt.inputs], modules=[muonScaleRes2017(),BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt", maxEntries=10000)
    if opt.year == "2017" and (opt.era in ["C"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017C(),BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017" and (opt.era in ["D"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017D(),BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017" and (opt.era in ["E"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017E(),BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017" and (opt.era in ["F"]):
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2017(),jmeCorrections_UL2017F(),BH2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018":
      p = PostProcessor(".", inputFiles(), modules=[muonScaleRes2018(),BH2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
  p.run()

if __name__ == "__main__":
    sys.exit(main())
