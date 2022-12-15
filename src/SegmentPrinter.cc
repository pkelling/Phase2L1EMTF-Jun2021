#include "L1Trigger/Phase2L1EMTF/interface/SegmentPrinter.h"

#include <algorithm>  // provides std::copy
#include <iostream>
#include <iterator>  // provides std::ostream_iterator
#include <vector>

#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"

#include "L1Trigger/Phase2L1EMTF/interface/Toolbox.h"

using namespace emtf::phase2;

struct SegmentPrinter::compress_sector {
  constexpr int operator()(int tp_endcap, int tp_sector) const {
    //tp_endcap = (tp_endcap == -1) ? 2 : tp_endcap;  // using endcap [1,2] convention
    //return ((tp_endcap - MIN_ENDCAP) * (NUM_TRIGSECTORS / 2)) + (tp_sector - MIN_TRIGSECTOR);
    return (tp_endcap * tp_sector); // -6..-1, 1..6
  }
};

struct SegmentPrinter::compress_cscid {
  constexpr int operator()(int tp_station, int tp_subsector, int tp_cscid, bool is_neighbor) const {
    if (not is_neighbor) {
      return (((tp_station == 1) ? (tp_subsector - 1) : 0) * 9) + (tp_cscid - 1);
    } else {
      return (tp_station == 1) ? (((tp_cscid - 1) / 3) + 18) : (((tp_cscid - 1) >= 3) + 9);
    }
  }
};

struct SegmentPrinter::pretty_print {
  std::ostream& operator()(std::ostream& os, const std::vector<int>& id_vec, const std::vector<int>& prop_vec) const {
    assert((not id_vec.empty()) and (not prop_vec.empty()));
    os << "(";
    std::copy(id_vec.begin(), --id_vec.end(), std::ostream_iterator<int>(os, ", "));
    std::copy(--id_vec.end(), id_vec.end(), std::ostream_iterator<int>(os, ")"));
    os << " ";
    std::copy(prop_vec.begin(), --prop_vec.end(), std::ostream_iterator<int>(os, " "));
    std::copy(--prop_vec.end(), prop_vec.end(), std::ostream_iterator<int>(os, ""));
    return os << std::endl;
  }
};

// _____________________________________________________________________________
void SegmentPrinter::print_impl(const csc_subsystem_tag::detid_type& detid,
                                const csc_subsystem_tag::digi_type& digi) const {
  static const int subsystem = L1TMuon::kCSC;
  static const int csc_bx_shift = -CSCConstants::LCT_CENTRAL_BX;

  int tp_subsector = toolbox::get_trigger_subsector(detid.station(), detid.chamber());

  std::vector<int> id_vec(7, 0);
  id_vec.at(0) = subsystem;
  id_vec.at(1) = detid.station();
  id_vec.at(2) = detid.ring();
  id_vec.at(3) = compress_sector{}(detid.endcap(), detid.triggerSector());
  id_vec.at(4) = compress_cscid{}(detid.station(), tp_subsector, digi.getCSCID(), false);
  id_vec.at(5) = 0;
  id_vec.at(6) = static_cast<int>(digi.getBX()) + csc_bx_shift;
  std::vector<int> prop_vec(4, 0);
  prop_vec.at(0) = digi.getStrip();
  prop_vec.at(1) = digi.getKeyWG();
  prop_vec.at(2) = digi.getPattern();
  prop_vec.at(3) = digi.isValid();
  pretty_print{}(std::cout, id_vec, prop_vec);

}

void SegmentPrinter::print_impl(const rpc_subsystem_tag::detid_type& detid,
                                const rpc_subsystem_tag::digi_type& digi) const {
  static const int subsystem = L1TMuon::kRPC;

  // Identifier for barrel RPC
  const bool is_barrel = (detid.region() == 0);

  // Identifier for iRPC (RE3/1, RE4/1)
  const bool is_irpc = ((not is_barrel) and (detid.station() >= 3) and (detid.ring() == 1));

  // Reject RPCb and RPCf-in-overlap-region (RE1/3, RE2/3)
  if (is_barrel or ((detid.station() <= 2) and (detid.ring() == 3)))
    return;

  int tp_chamber = ((detid.sector() - 1) * (is_irpc ? 3 : 6)) + detid.subsector();
  int tp_sector = toolbox::get_trigger_sector(detid.ring(), detid.station(), tp_chamber);
  int tp_cscid = toolbox::get_trigger_cscid(detid.ring(), detid.station(), tp_chamber);
  int tp_subsector = toolbox::get_trigger_subsector(detid.station(), tp_chamber);

  std::vector<int> id_vec(7, 0);
  id_vec.at(0) = subsystem;
  id_vec.at(1) = detid.station();
  id_vec.at(2) = detid.ring();
  id_vec.at(3) = compress_sector{}(detid.region(), tp_sector);
  id_vec.at(4) = compress_cscid{}(detid.station(), tp_subsector, tp_cscid, false);
  id_vec.at(5) = detid.layer();
  id_vec.at(6) = digi.BunchX();
  std::vector<int> prop_vec(4, 0);
  prop_vec.at(0) = digi.firstClusterStrip();
  prop_vec.at(1) = detid.roll();
  prop_vec.at(2) = digi.clusterSize();
  prop_vec.at(3) = true;
  pretty_print{}(std::cout, id_vec, prop_vec);
}

void SegmentPrinter::print_impl(const gem_subsystem_tag::detid_type& detid,
                                const gem_subsystem_tag::digi_type& digi) const {
  static const int subsystem = L1TMuon::kGEM;

  int tp_sector = toolbox::get_trigger_sector(detid.ring(), detid.station(), detid.chamber());
  int tp_cscid = toolbox::get_trigger_cscid(detid.ring(), detid.station(), detid.chamber());
  int tp_subsector = toolbox::get_trigger_subsector(detid.station(), detid.chamber());

  std::vector<int> id_vec(7, 0);
  id_vec.at(0) = subsystem;
  id_vec.at(1) = detid.station();
  id_vec.at(2) = detid.ring();
  id_vec.at(3) = compress_sector{}(detid.region(), tp_sector);
  id_vec.at(4) = compress_cscid{}(detid.station(), tp_subsector, tp_cscid, false);
  id_vec.at(5) = detid.layer();
  id_vec.at(6) = digi.bx();
  std::vector<int> prop_vec(4, 0);
  prop_vec.at(0) = digi.pads().front();
  prop_vec.at(1) = detid.roll();
  prop_vec.at(2) = digi.pads().back() - digi.pads().front() + 1;
  prop_vec.at(3) = digi.isValid();
  pretty_print{}(std::cout, id_vec, prop_vec);

  std::cout << "GE" << detid.station() << "1-" << id_vec.at(4) << "\t" << digi.pads().front() << ", " << detid.roll() << ", " << prop_vec.at(2) << std::endl;

}

void SegmentPrinter::print_impl(const me0_subsystem_tag::detid_type& detid,
                                const me0_subsystem_tag::digi_type& digi) const {
  static const int subsystem = L1TMuon::kME0;
  static const int me0_bx_shift = -CSCConstants::LCT_CENTRAL_BX;

  int fake_tp_station = 2;  // ME0 is 20-deg chamber like ME2/1
  int fake_tp_ring = 1;
  int tp_sector = toolbox::get_trigger_sector(fake_tp_ring, fake_tp_station, detid.chamber());
  int tp_cscid = toolbox::get_trigger_cscid(fake_tp_ring, fake_tp_station, detid.chamber());
  int tp_subsector = toolbox::get_trigger_subsector(fake_tp_station, detid.chamber());

  std::vector<int> id_vec(7, 0);
  id_vec.at(0) = subsystem;
  id_vec.at(1) = detid.station();
  id_vec.at(2) = 4;
  id_vec.at(3) = compress_sector{}(detid.region(), tp_sector);
  id_vec.at(4) = compress_cscid{}(fake_tp_station, tp_subsector, tp_cscid, false);
  id_vec.at(5) = detid.layer();
  id_vec.at(6) = static_cast<int>(digi.getBX()) + me0_bx_shift;
  std::vector<int> prop_vec(4, 0);
  prop_vec.at(0) = digi.getPhiposition();
  prop_vec.at(1) = digi.getPartition();
  prop_vec.at(2) = static_cast<int>(digi.getDeltaphi()) * (static_cast<int>(digi.getBend()) * 2 - 1);
  prop_vec.at(3) = digi.isValid();
  pretty_print{}(std::cout, id_vec, prop_vec);

  std::cout << "ME0-" << detid.chamber() << "\t" << digi.getStrip() << ", " << digi.getPartition() << ", " << std::endl;

}

void SegmentPrinter::print_impl(const EMTFHit& hit) const {
  std::vector<int> id_vec(7, 0);
  id_vec.at(0) = hit.subsystem();
  id_vec.at(1) = hit.station();
  id_vec.at(2) = hit.ring();
  id_vec.at(3) = compress_sector{}(hit.endcap(), hit.sector());
  id_vec.at(4) = compress_cscid{}(hit.station(), hit.subsector(), hit.cscid(), hit.neighbor());
  id_vec.at(5) = hit.gemdl();
  id_vec.at(6) = hit.bx();
  std::vector<int> prop_vec(15, 0);
  prop_vec.at(0) = hit.emtfChamber();
  prop_vec.at(1) = hit.emtfSegment();
  prop_vec.at(2) = hit.emtfPhi();
  prop_vec.at(3) = hit.emtfBend();
  prop_vec.at(4) = hit.emtfTheta1();
  prop_vec.at(5) = hit.emtfTheta2();
  prop_vec.at(6) = hit.emtfQual1();
  prop_vec.at(7) = hit.emtfQual2();
  prop_vec.at(8) = hit.emtfTime();
  prop_vec.at(9) = hit.zones();
  prop_vec.at(10) = hit.timezones();
  prop_vec.at(11) = hit.cscfr();
  prop_vec.at(12) = hit.gemdl();
  prop_vec.at(13) = hit.bx();
  prop_vec.at(14) = hit.valid();
  pretty_print{}(std::cout, id_vec, prop_vec);
}

void SegmentPrinter::print_impl(const EMTFTrack& trk) const {
  std::vector<int> id_vec(3, 0);
  id_vec.at(0) = compress_sector{}(trk.endcap(), trk.sector());
  id_vec.at(1) = trk.bx();
  id_vec.at(2) = trk.unconstrained();
  std::vector<int> prop_vec(4, 0);
  prop_vec.at(0) = trk.modelInvpt();
  prop_vec.at(1) = trk.modelPhi();
  prop_vec.at(2) = trk.modelEta();
  prop_vec.at(3) = trk.modelQual();
  pretty_print{}(std::cout, id_vec, prop_vec);
}
