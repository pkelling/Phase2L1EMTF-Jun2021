#include "L1Trigger/Phase2L1EMTF/interface/SectorProcessor.h"

#include <algorithm>  // provides std::find
#include <array>
#include <iostream>
#include <iterator>  // provides std::make_move_iterator
#include <type_traits>
#include <map>
#include <utility>
#include <variant>
#include <vector>

#include "L1Trigger/Phase2L1EMTF/interface/EMTFWorker.h"
#include "L1Trigger/Phase2L1EMTF/interface/EMTFModel.h"
#include "L1Trigger/Phase2L1EMTF/interface/GeometryHelper.h"
#include "L1Trigger/Phase2L1EMTF/interface/ConditionHelper.h"
#include "L1Trigger/Phase2L1EMTF/interface/SegmentFormatter.h"
#include "L1Trigger/Phase2L1EMTF/interface/SegmentPrinter.h"
#include "L1Trigger/Phase2L1EMTF/interface/TrackFormatter.h"
#include "L1Trigger/Phase2L1EMTF/interface/Toolbox.h"

//#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

using namespace emtf::phase2;

void SectorProcessor::process(const EMTFWorker& iWorker,
                              int endcap,
                              int sector,
                              const edm::EventID& evt_id,
                              const SubsystemCollection& muon_primitives,
                              EMTFHitCollection& out_hits,
                              EMTFTrackCollection& out_tracks) const {

  // Loop over BX
  for (int bx = iWorker.minBX_; bx <= iWorker.maxBX_; ++bx) {
    // 1 - Preprocessing
    EMTFHitCollection sector_hits;
    process_step_1(iWorker, endcap, sector, bx, muon_primitives, sector_hits);

    // 2 - Real processing
    // Only BX=0 is supported at the moment
    EMTFTrackCollection sector_tracks;
    if (bx == 0) {
      process_step_2(iWorker, endcap, sector, bx, sector_hits, sector_tracks);
    }

    // 3 - Postprocessing
    out_hits.insert(
        out_hits.end(), std::make_move_iterator(sector_hits.begin()), std::make_move_iterator(sector_hits.end()));
    out_tracks.insert(
        out_tracks.end(), std::make_move_iterator(sector_tracks.begin()), std::make_move_iterator(sector_tracks.end()));

#define EMTF_DUMP_INFO
#ifdef EMTF_DUMP_INFO
    // Dump debugging info after the last sector only
    if ((endcap == MAX_ENDCAP) and (sector == MAX_TRIGSECTOR) and (bx == 0)) {
      dump_input_output(evt_id, muon_primitives, out_hits, out_tracks);
    }
#endif  // EMTF_DUMP_INFO is defined

  }  // end loop over BX

}

void SectorProcessor::process_step_1(const EMTFWorker& iWorker,
                                     int endcap,
                                     int sector,
                                     int bx,
                                     const SubsystemCollection& muon_primitives,
                                     EMTFHitCollection& sector_hits) const {
  // For CSC, keep a list of wire ambiguity. Store the list in a map with key: (detid, bx), value: (wire,).
  // For GEM, keep a list of coincidence pads. Store the list in a map with key: (detid, bx), value: (roll, pad_lo, pad_hi).
  // For RPC and ME0, do nothing.
  std::map<std::pair<uint32_t, uint16_t>, std::vector<uint16_t> > csc_chamber_wire_ambi;
  std::map<std::pair<uint32_t, uint16_t>, std::vector<std::array<uint16_t, 3> > > gem_chamber_copad_vec;

  // Loop over muon_primitives (1st pass)
  for (const auto& [a, b, c] : muon_primitives) {
    // clang-format off
    std::visit([&](auto&& subsystem, auto&& detid, auto&& digi) {
      using T1 = std::decay_t<decltype(subsystem)>;
      using T2 = std::decay_t<decltype(detid)>;
      using T3 = std::decay_t<decltype(digi)>;
      // Enable if (subsystem, detid, digi) are consistent
      if constexpr (std::is_same_v<typename T1::detid_type, T2> and std::is_same_v<typename T1::digi_type, T3>) {
        // Enable if subsystem is CSC/RPC/GEM/ME0
        if constexpr (std::is_same_v<T1, csc_subsystem_tag>) {
          uint16_t tp_wire = digi.getKeyWG();
          auto akey = std::make_pair(detid.rawId(), digi.getBX());
          // If key and value both exist, do nothing. If key exists, but not value, insert value.
          // If neither key nor value exists, insert both key and value.
          auto found = csc_chamber_wire_ambi.find(akey);
          if (found != csc_chamber_wire_ambi.end()) {
            auto inner_found = std::find(found->second.begin(), found->second.end(), tp_wire);
            if (inner_found != found->second.end()) {
              // Do nothing
            } else {
              found->second.push_back(tp_wire);
            }
          } else {
            csc_chamber_wire_ambi[akey].push_back(tp_wire);
          }

        } else if constexpr (std::is_same_v<T1, rpc_subsystem_tag>) {
          // Do nothing

        } else if constexpr (std::is_same_v<T1, gem_subsystem_tag>) {
          uint16_t tp_layer = detid.layer();
          uint16_t tp_roll = detid.roll();
          uint16_t tp_pad_lo = digi.pads().front();
          uint16_t tp_pad_hi = digi.pads().back();
          bool tp_valid = digi.isValid();
          // Remove layer number and roll number from detid
          gem_subsystem_tag::detid_type detid_mod(
              detid.region(), detid.ring(), detid.station(), 0, detid.chamber(), 0);
          auto akey = std::make_pair(detid_mod.rawId(), digi.bx());
          if (tp_valid and (tp_layer == 1)) {  // layer 1 is used as incidence
            // If key does not exist, insert an empty vector. If key exists, do nothing.
            decltype(gem_chamber_copad_vec)::mapped_type avec;
            gem_chamber_copad_vec.insert({akey, avec});
          } else if (tp_valid and (tp_layer == 2)) {  // layer 2 is used as coincidence
            // If key does not exist, insert an empty vector and push a value. If key exists, push a value.
            decltype(gem_chamber_copad_vec)::mapped_type::value_type aval{{tp_roll, tp_pad_lo, tp_pad_hi}};
            gem_chamber_copad_vec[akey].push_back(std::move(aval));
          }

        } else if constexpr (std::is_same_v<T1, me0_subsystem_tag>) {
          // Do nothing
        } else {
          // Make sure every subsystem type has been visited
          static_assert(dependent_false<T1>::value, "unreachable!");
        }         // end inner constexpr if statement
      }           // end outer constexpr if statement
    }, a, b, c);  // end visit
    // clang-format on

  }  // end loop

  // Convert/format input segments
  SegmentFormatter formatter;
  EMTFHitCollection substitutes;

  // Loop over muon_primitives (2nd pass)
  for (const auto& [a, b, c] : muon_primitives) {
    SegmentFormatter::ChamberInfo chminfo;
    EMTFHit hit;
    int strategy = 0;  // default strategy

    // clang-format off
    std::visit([&](auto&& subsystem, auto&& detid, auto&& digi) {
      using T1 = std::decay_t<decltype(subsystem)>;
      using T2 = std::decay_t<decltype(detid)>;
      using T3 = std::decay_t<decltype(digi)>;
      using T4 = typename T1::detgeom_type;
      // Enable if (subsystem, detid, digi) are consistent
      if constexpr (std::is_same_v<typename T1::detid_type, T2> and std::is_same_v<typename T1::digi_type, T3>) {
        // Enable if subsystem is CSC/RPC/GEM/ME0
        if constexpr (std::is_same_v<T1, csc_subsystem_tag>) {
          // For CSC, send the list of wire ambiguity
          auto akey = std::make_pair(detid.rawId(), digi.getBX());
          chminfo.wire_ambi = csc_chamber_wire_ambi.at(akey);

        } else if constexpr (std::is_same_v<T1, rpc_subsystem_tag>) {
          // Do nothing

        } else if constexpr (std::is_same_v<T1, gem_subsystem_tag>) {
          // For GEM, send the list of coincidence pads
          // Remove layer number and roll number from detid
          gem_subsystem_tag::detid_type detid_mod(
              detid.region(), detid.ring(), detid.station(), 0, detid.chamber(), 0);
          auto akey = std::make_pair(detid_mod.rawId(), digi.bx());
          chminfo.copad_vec = gem_chamber_copad_vec.at(akey);

        } else if constexpr (std::is_same_v<T1, me0_subsystem_tag>) {
          // Do nothing

        } else {
          // Make sure every subsystem type has been visited
          static_assert(dependent_false<T1>::value, "unreachable!");
        }         // end inner constexpr if statement

        // Get subsystem geometry and do the conversion
        auto&& detgeom = iWorker.geom_helper_->get<T4>();

        formatter.format(endcap, sector, bx, strategy, detgeom, detid, digi, chminfo, hit);

        // Try again with a different strategy
        if (not hit.valid()) {
          strategy++;
          formatter.format(endcap, sector, bx, strategy, detgeom, detid, digi, chminfo, hit);
        }
      }           // end outer constexpr if statement
    }, a, b, c);  // end visit
    // clang-format on

    // Does not belong to this sector
    if (not hit.valid())
      continue;

    // Keep the valid segment
    if (strategy == 0) {
      sector_hits.push_back(std::move(hit));
    } else {
      substitutes.push_back(std::move(hit));
    }
  }  // end loop

  // Insert substitutes at the end of sector_hits
  sector_hits.insert(
      sector_hits.end(), std::make_move_iterator(substitutes.begin()), std::make_move_iterator(substitutes.end()));

  // Modify the segments in-place
  std::map<int, int> counter_map;  // key: emtf_chamber

  for (auto&& hit : sector_hits) {
    // Count num of valid segments in this chamber
    int emtf_chamber = hit.emtfChamber();
    counter_map[emtf_chamber]++;                          // increment counter
    int emtf_segment = counter_map.at(emtf_chamber) - 1;  // decrement to get the original count
    hit.setEmtfSegment(emtf_segment);                     // assign emtf_segment
  }
}

void SectorProcessor::process_step_2(const EMTFWorker& iWorker,
                                     int endcap,
                                     int sector,
                                     int bx,
                                     const EMTFHitCollection& sector_hits,
                                     EMTFTrackCollection& sector_tracks) const {
  // Exit early if sector is empty
  bool early_exit = sector_hits.empty();

  if (early_exit)
    return;

  const NdArrayDesc& input_shape = iWorker.model_->get_input_shape();
  const NdArrayDesc& output_shape = iWorker.model_->get_output_shape();
  const unsigned num_segments = iWorker.model_->get_num_segments();
  const unsigned num_tracks = iWorker.model_->get_num_tracks();

  // Model input and output
  std::vector<int> in0(input_shape.num_elements(), 0);
  std::vector<int> out(output_shape.num_elements(), 0);

  // Fill values
  for (auto&& hit : sector_hits) {
    const int emtf_chamber = hit.emtfChamber();
    const int emtf_segment = hit.emtfSegment();
    emtf_assert(hit.valid() == true);  // segment must be valid

    // Accept at most 2 segments
    if (not(static_cast<unsigned>(emtf_segment) < num_segments))
      continue;

    // Populate the variables
    //
    // +-------------+-------------+-------------+-------------+
    // | emtf_phi    | emtf_bend   | emtf_theta1 | emtf_theta2 |
    // +-------------+-------------+-------------+-------------+
    // | emtf_qual1  | emtf_qual2  | emtf_time   | seg_zones   |
    // +-------------+-------------+-------------+-------------+
    // | seg_tzones  | seg_cscfr   | seg_gemdl   | seg_bx      |
    // +-------------+-------------+-------------+-------------+
    // | seg_valid   |             |             |             |
    // +-------------+-------------+-------------+-------------+

    const unsigned iseg = (emtf_chamber * num_segments) + emtf_segment;
    const unsigned ivar = 0;
    auto in0_iter = std::next(in0.begin(), input_shape.get_index({iseg, ivar}));

    *(in0_iter++) = hit.emtfPhi();
    *(in0_iter++) = hit.emtfBend();
    *(in0_iter++) = hit.emtfTheta1();
    *(in0_iter++) = hit.emtfTheta2();
    *(in0_iter++) = hit.emtfQual1();
    *(in0_iter++) = hit.emtfQual2();
    *(in0_iter++) = hit.emtfTime();
    *(in0_iter++) = hit.zones();
    *(in0_iter++) = hit.timezones();
    *(in0_iter++) = hit.cscfr();
    *(in0_iter++) = hit.gemdl();
    *(in0_iter++) = hit.bx();
    *(in0_iter++) = hit.valid();
  }  // end loop

  // Fit
  iWorker.model_->fit(in0, out);

  // Convert/format output tracks
  TrackFormatter formatter;
  const unsigned model_version = iWorker.model_->version();
  const bool unconstrained = iWorker.model_->unconstrained();

  // Extract results
  for (unsigned itrk = 0, ivar = 0; itrk < num_tracks; ++itrk) {
    EMTFTrack trk;

    // Get the span of data and do the conversion
    auto out_iter = std::next(out.begin(), output_shape.get_index({itrk, ivar}));
    auto out_iter_end =
        ((itrk + 1) < num_tracks) ? std::next(out.begin(), output_shape.get_index({itrk + 1, ivar})) : out.end();
    std::vector<int> trk_data(out_iter, out_iter_end);
    formatter.format(endcap, sector, bx, model_version, unconstrained, trk_data, trk);

    // Skip the invalid track
    if (not trk.valid())
      continue;

    // Keep the valid track
    sector_tracks.push_back(std::move(trk));
  }  // end loop
}

void SectorProcessor::dump_input_output(const edm::EventID& evt_id,
                                        const SubsystemCollection& muon_primitives,
                                        const EMTFHitCollection& out_hits,
                                        const EMTFTrackCollection& out_tracks) const {
  const std::string bold_seq = "\033[1m";
  const std::string reset_seq = "\033[0m";

  std::cout << "Processing " << evt_id << std::endl;
  std::cout << "[" << bold_seq << "RX" << reset_seq << "]" << std::endl;

  SegmentPrinter printer;

  // Loop over muon_primitives
  for (const auto& [a, b, c] : muon_primitives) {
    // clang-format off
    std::visit([&](auto&& subsystem, auto&& detid, auto&& digi) {
      using T1 = std::decay_t<decltype(subsystem)>;
      using T2 = std::decay_t<decltype(detid)>;
      using T3 = std::decay_t<decltype(digi)>;
      // Enable if (subsystem, detid, digi) are consistent
      if constexpr (std::is_same_v<typename T1::detid_type, T2> and std::is_same_v<typename T1::digi_type, T3>) {
        printer.print(detid, digi);
      }           // end constexpr if statement
    }, a, b, c);  // end visit
    // clang-format on

  }  // end loop

  std::cout << "[" << bold_seq << "TX#0" << reset_seq << "]" << std::endl;

  // Loop over converted EMTF hits
  for (auto&& hit : out_hits) {
    printer.print(hit);
  }  // end loop

  std::cout << "[" << bold_seq << "TX#1" << reset_seq << "]" << std::endl;

  // Loop over EMTF tracks
  for (auto&& trk : out_tracks) {
    printer.print(trk);
  }  // end loop
}


