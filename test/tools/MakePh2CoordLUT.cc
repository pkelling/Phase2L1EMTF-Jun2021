#include <cmath>
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


class MakePh2CoordLUT : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit MakePh2CoordLUT(const edm::ParameterSet&);
  ~MakePh2CoordLUT() override;

private:

  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void endRun(const edm::Run&, const edm::EventSetup&) override;

  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  // Generate LUTs
  void generateLUTs();

private:
  const edm::ParameterSet config_;

  int verbose_;

  bool done_;

};

// _____________________________________________________________________________

MakePh2CoordLUT::MakePh2CoordLUT(const edm::ParameterSet& iConfig)
    : config_(iConfig),
      done_(false){
  std::cout << "Running MakeCoordLUT\n";
}

MakePh2CoordLUT::~MakePh2CoordLUT() {}

void MakePh2CoordLUT::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  /// Geometry setup
  std::cout << "Beginning Run\n";
}

void MakePh2CoordLUT::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void MakePh2CoordLUT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (done_)
    return;

  generateLUTs();

  done_ = true;
  return;
}

// _____________________________________________________________________________
void MakePh2CoordLUT::generateLUTs() {
  std::cout << "Running generateLUTs()\n";
}

// DEFINE THIS AS A PLUG-IN
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakePh2CoordLUT);
