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


#define EMTF_CREATE_FW_LUTS
#ifdef EMTF_CREATE_FW_LUTS
  if( endcap == 1 && sector == 1){ // only runs once
    //create_me0_luts(iWorker);
    //create_gem_luts(iWorker);
    //test_dump(iWorker, muon_primitives);
    create_csc_luts(iWorker);
  }
  return;
#endif

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

//#define EMTF_DUMP_INFO
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



void SectorProcessor::test_dump(const EMTFWorker& iWorker, const SubsystemCollection& muon_primitives) const {
  
  // Loop over muon_primitives
  for (const auto& [a, b, c] : muon_primitives) {
    // clang-format off
    std::visit([&](auto&& subsystem, auto&& detid, auto&& digi) {
      using T1 = std::decay_t<decltype(subsystem)>;
      using T2 = std::decay_t<decltype(detid)>;
      using T3 = std::decay_t<decltype(digi)>;
      // Enable if (subsystem, detid, digi) are consistent
      if constexpr (std::is_same_v<typename T1::detid_type, T2> and std::is_same_v<typename T1::digi_type, T3>) {
        if constexpr (std::is_same_v<T1, gem_subsystem_tag>) {
          SegmentFormatter::ChamberInfo chminfo;
          EMTFHit hit_o;

          // Coincident layer features
          chminfo.copad_vec = { {{8, 26, 28}}, {{7,25,25}} }; // { {roll, pad_start, pad_stop} }

          // digi features
          GEMPadDigiCluster digiTest(digi.pads(), digi.bx());
          std::cout << digi << std::endl;

          // detid features
          int chamber = 3;
          int roll = 8;
          //(region, ring, station, layer, chamber, roll)
          GEMDetId detidTest(1, 1, 2, 1, chamber, roll);


          using T4 = typename T1::detgeom_type;
          auto&& detgeom = iWorker.geom_helper_->get<T4>();

          SegmentFormatter formatter;
                        //ecap, sect
          formatter.format(1, 1, 0, 0, detgeom, detidTest, digiTest, chminfo, hit_o);
          std::cout << "Strip " << hit_o.strip() << " - Phi: " << hit_o.emtfPhi() << std::endl;
          std::cout << "WG  " << hit_o.wire1() << " - Theta: " << hit_o.emtfTheta1() << std::endl << std::endl;

        }
      }           // end constexpr if statement
    }, a, b, c);  // end visit
    // clang-format on

  }  // end loop

  
}

int SectorProcessor::ge11_phi_conv(int chamber, int ph_init, int hs) const {

  bool  ph_reverse = (chamber % 2 == 1) ? 1 : 0;

  int factor = 3256 / 2;

  // set ph
  int mult = hs * factor;
  int ph_tmp = mult >> 10;
  int ph; 
  if (ph_reverse) ph = ph_init - ph_tmp;
  else            ph = ph_init + ph_tmp;
       
  std::cout << hs << "\t" << ph << std::endl; // 3256
  return ph;
}


int SectorProcessor::ge21_phi_conv(int chamber, int ph_init, int hs) const {

  bool  ph_reverse = (chamber % 2 == 1) ? 1 : 0;

  int factor = 1672; // 1669; // 3337.2 / 2 

  // set ph
  int mult = hs * factor;
  int ph_tmp = mult >> 10;
  int ph; 
  //if (ph_reverse) ph = ph_init - ph_tmp;
  //else            ph = ph_init + ph_tmp;
       
  if (ph_reverse){ 
    ph = ph_init - ph_tmp;
    /*
    switch(hs){
      case 0 ... 26:
        ph-=2; break;
      case 27 ... 57:
        ph-=1; break;
        break;
      case 177 ... 266:
        ph+=1; break;
      case 485 ... 636:
        ph-=1; break;
      case 727 ... 756:
        ph+=1; break;
      case 757 ... 766:
        ph+=2; break;
    }
    */
  }
  else{            
    ph = ph_init + ph_tmp;
    /*
    switch(hs){
      case 0 ... 26:
        ph+=2; break;
      case 27 ... 57:
        ph+=1; break;
        break;
      case 177 ... 266:
        ph-=1; break;
      case 485 ... 636:
        ph+=1; break;
      case 727 ... 756:
        ph-=1; break;
      case 757 ... 766:
        ph-=2; break;
    }
    */
  }

  //std::cout << hs << "\t" << ph << std::endl;
  return ph;
}


void SectorProcessor::create_gem_luts(const EMTFWorker& iWorker) const {
  
  SegmentFormatter::ChamberInfo chminfo;
  SegmentFormatter formatter;
  EMTFHit hit_o;
  
  // detid features
  int region = 1; // +- 1 for endcap
  uint16_t roll = 7;
  int chamber_start, chamber_stop, max_pad;

  int station = 2; // ge11 or ge21

  if(station == 1){ //ge11
    chamber_start = 2; //neighbor chamber of sector 1
    chamber_stop = 8; // 8
    max_pad = 191;
  }
  else{ // ge21
    chamber_start = 1; //neighbor chamber of sector 1
    chamber_stop = 4; //4
    max_pad = 383;
  }


  // Test SW vs FW implementations
  int ph_init;
  int ph;
  int halfpad, cl_sz;
  int max_error = 0;
  int total_error = 0;
  int count = 0;

  for(int chamber = chamber_start; chamber <= chamber_stop; chamber++){
    GEMDetId detidTest(region, 1, station, 1, chamber, roll); // (region, ring, station, layer, chamber, roll )

    using T4 = typename gem_subsystem_tag::detgeom_type;
    auto&& detgeom = iWorker.geom_helper_->get<T4>();

    std::cout << "\nChamber " << chamber << std::endl;
    for(uint16_t i=0; i<=max_pad; i++){
      
      chminfo.copad_vec = { {{roll, i, i}} }; // {{roll, pad_start, pad_stop}}

      // Set incident pads (layer 1) and coincident pads (layer 2)
      cl_sz = 0;
      halfpad = (i << 1) + cl_sz;

      std::vector<uint16_t> digi_pads = {i};
      GEMPadDigiCluster digiTest(digi_pads, 0); // (digi_pads, bx)

      formatter.format(1, 1, 0, 0, detgeom, detidTest, digiTest, chminfo, hit_o); // (endcap, sector, ...)
      //std::cout << (hit_o.strip()<<1) << "\t" << hit_o.emtfPhi() << std::endl;
      //std::cout << hit_o.emtfPhi() << "," << std::endl;

      //if(i == 0)
        //ph_init = hit_o.emtfPhi();

      if(station == 1)
        ph = ge11_phi_conv(chamber, ph_init, halfpad);
      else{
        switch(chamber){
          case 1:
            ph_init = 1345;
            break;
          case 2:
            ph_init = 1295;
            break;
          case 3:
            ph_init = 3745;
            break;
          case 4:
            ph_init = 3695;
            break;
        }
        ph = ge21_phi_conv(chamber, ph_init, halfpad); 
      }

      total_error += std::abs(ph - hit_o.emtfPhi());
      if( std::abs(ph - hit_o.emtfPhi()) > max_error)
        max_error = std::abs(ph - hit_o.emtfPhi());
      if( std::abs(ph - hit_o.emtfPhi()) >= 2){
        std::cout << "Pad: " << i << "\t" << (ph - hit_o.emtfPhi()) << std::endl; 
        //std::cout << (ph - hit_o.emtfPhi()) << std::endl;
        count++;
      }


      // Testing Even number of pads
      if(i<(max_pad)){
        cl_sz = 1; // 6 pads
        halfpad = (i << 1) + cl_sz;

        std::vector<uint16_t> digi_pads2 = {i,uint16_t(i+1)};
        GEMPadDigiCluster digiTest2(digi_pads2, 0); // (digi_pads, bx)

        formatter.format(1, 1, 0, 0, detgeom, detidTest, digiTest2, chminfo, hit_o); // (endcap, sector, ...)
        //std::cout << ((hit_o.strip()<<1)+1) << "\t" << hit_o.emtfPhi() << std::endl;
        //std::cout << hit_o.emtfPhi() << "," <<  std::endl;

        if(station == 1)
          ph = ge11_phi_conv(chamber, ph_init, halfpad);
        else{
          switch(chamber){
            case 1:
              ph_init = 1345;
              break;
            case 2:
              ph_init = 1295;
              break;
            case 3:
              ph_init = 3745;
              break;
            case 4:
              ph_init = 3695;
              break;
          }
          ph = ge21_phi_conv(chamber, ph_init, halfpad); 
        }
        
        total_error += std::abs(ph - hit_o.emtfPhi());
        if( std::abs(ph - hit_o.emtfPhi()) > max_error)
          max_error = std::abs(ph - hit_o.emtfPhi());
        if( std::abs(ph - hit_o.emtfPhi()) >= 2){
          std::cout << "Pad mid: " << i << "\t" << (ph - hit_o.emtfPhi()) << std::endl;
          //std::cout << (ph - hit_o.emtfPhi()) << std::endl;
          count++;
        }
      }
    }
  }
  std::cout << "Max Error: " << max_error << std::endl;
  std::cout << "Total Error: " << total_error << std::endl;
  std::cout << "Count Bad: " << count << std::endl;

  /*
  // Print ph_init values
  for(int chamber = chamber_start; chamber <= chamber_stop; chamber++){
    roll = 7;
    GEMDetId detidTest(region, 1, station, 1, chamber, roll); // (region, ring, station, layer, chamber, roll )

    using T4 = typename gem_subsystem_tag::detgeom_type;
    auto&& detgeom = iWorker.geom_helper_->get<T4>();

    uint16_t i=max_pad; //quick compiler error fix
    i=0;

    chminfo.copad_vec = { {{roll, i, i}} }; // {{roll, pad_start, pad_stop}}
    std::vector<uint16_t> digi_pads = {i};
    GEMPadDigiCluster digiTest(digi_pads, 0); // (digi_pads, bx)

    formatter.format(1, 1, 0, 0, detgeom, detidTest, digiTest, chminfo, hit_o); // (endcap, sector, ...)
    std::cout << "Chamber: " << chamber << std::hex << "\t" << hit_o.emtfPhi() << std::dec << std::endl;
  }
  */


  /*
  // print theta LUTs for sector
  for(int chamber = chamber_start; chamber <= chamber_stop; chamber++){
    std::cout << "\nChamber " << chamber << std::endl;
      
    for(roll=1; roll<=8; roll++){
      GEMDetId detidTest(region, 1, station, 1, chamber, roll); // (region, ring, station, layer, chamber, roll )

      using T4 = typename gem_subsystem_tag::detgeom_type;
      auto&& detgeom = iWorker.geom_helper_->get<T4>();

      uint16_t i=max_pad;
      chminfo.copad_vec = { {{roll, i, i}} }; // {{roll, pad_start, pad_stop}}
      std::vector<uint16_t> digi_pads = {i};
      GEMPadDigiCluster digiTest(digi_pads, 0); // (digi_pads, bx)

      formatter.format(1, 1, 0, 0, detgeom, detidTest, digiTest, chminfo, hit_o); // (endcap, sector, ...)
      std::cout << std::hex << hit_o.emtfTheta1() << std::endl;
    }
    std::cout << std::endl;
  }
  */

}


void SectorProcessor::create_csc_luts(const EMTFWorker& iWorker) const {

  SegmentFormatter formatter;
  int endcap = 1; // +1
  int station = 1;
  int ring = 1;
  int layer = 0;
  int chamber = 27;
  uint16_t strip = 15;
  uint16_t quality = 14;
  uint16_t keywire = 10;
  uint16_t pattern = 8;
  uint16_t bend = 1;
  uint16_t bx = 8; // for DIGI

  //unsigned low_point_phi = 0;
  
  for( chamber = 1; chamber < 36; chamber++){
    for( strip=0; strip < 1; strip++){
      for(keywire = 0; keywire<11; keywire+=10){

        EMTFHit hit_o;
        CSCDetId detid(endcap, station, ring, chamber, layer);
        int sector_i = toolbox::get_trigger_sector(detid.ring(), detid.station(), detid.chamber());

        CSCCorrelatedLCTDigi digi(1, 1, quality, keywire, strip, pattern, bend, bx, 1);
        digi.setType(1);

        uint16_t valid_cscid = toolbox::get_trigger_cscid(ring, station, chamber);
        digi.setCSCID(valid_cscid);

        using T4 = typename csc_subsystem_tag::detgeom_type;
        auto&& detgeom = iWorker.geom_helper_->get<T4>();
        SegmentFormatter::ChamberInfo chminfo;

        std::map<std::pair<uint32_t, uint16_t>, std::vector<uint16_t> > csc_chamber_wire_ambi;
        auto akey = std::make_pair(detid.rawId(), digi.getBX());
        uint16_t tp_wire = digi.getKeyWG();
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
        chminfo.wire_ambi = csc_chamber_wire_ambi.at(akey);

        formatter.format(endcap, sector_i, 0, 0, detgeom, detid, digi, chminfo, hit_o);

        std::cout << detid.chamberName() << " - " << hit_o.strip() << "," << hit_o.wire1() << "\tPhi: " << hit_o.emtfPhi() << std::endl;
        
        //std::cout << "St: " << station << "\tRing: " << ring << "\tSector: " << sector_i << std::endl;
        //std::cout << "Strip " << hit_o.strip() << " - Phi: " << hit_o.emtfPhi() << std::endl;
        //std::cout << "WG  " << hit_o.wire1() << " - Theta: " << hit_o.emtfTheta1() << std::endl << std::endl;

      }
      std::cout << std::endl;
    }
  }

}


void SectorProcessor::create_me0_luts(const EMTFWorker& iWorker) const {
  static const int me0_max_partition = 8;  // limited to eta of 2.4
  //static const int me0_nstrips = 384;
  //static const int me0_nphipositions = me0_nstrips * 2;

  // create detid, detgeom, and digi
  int region = +1; // endcap +- 1
  int chamber = 1;
  int roll = 1;
  int hs = 75;

  std::cout << std::endl;
  for(chamber = 1; chamber<=15; chamber++){
  for(roll = 0; roll < me0_max_partition; roll++){
    std::cout << "Chamber: " << chamber << "\tHalfstrip: " << hs << std::endl;
    SegmentFormatter::ChamberInfo chminfo;
    SegmentFormatter formatter;
    EMTFHit hit_o;
    
    ME0DetId detid(region, 0, chamber, roll); // (region, layer, chamber, roll )

    using T4 = typename me0_subsystem_tag::detgeom_type;
    auto&& detgeom = iWorker.geom_helper_->get<T4>();

    //std::vector<uint16_t> digi_pads = {i};
    //GEMPadDigiCluster digiTest(digi_pads, 0); // (digi_pads, bx)
    ME0TriggerDigi digi(chamber, 5, hs, roll, 20, 1, 8);
    /*ME0TriggerDigi digi(const int ichamberid,
                                 const int iquality,
                                 const int iphiposition,
                                 const int ipartition,
                                 const int ideltaphi,
                                 const int ibend,
                                 const int ibx)
    */

    int sector;
    if( chamber == 1)
      sector = 6;
    else if(chamber < 5)
      sector = 1;
    else if(chamber < 8)
      sector = 2;
    else if(chamber < 11)
      sector = 3;
    else if(chamber < 14)
      sector = 4;
    else
      sector = 5;

    formatter.format(1,sector,0,0,detgeom, detid, digi, chminfo, hit_o);
    /*SegmentFormatter::format_impl(int endcap, 1:+  2:-
                                     int sector,
                                     int bx,
                                     int strategy,
                                     const ME0Geometry& detgeom,
                                     const me0_subsystem_tag::detid_type& detid,
                                     const me0_subsystem_tag::digi_type& digi,
                                     const ChamberInfo& chminfo,
                                     EMTFHit& hit) const {
    */
    //std::cout << "Strip " << hit_o.strip() << " - Phi: " << hit_o.emtfPhi() << std::endl << std::endl;
    //std::cout << "WG  " << hit_o.wire1() << " - Theta: " << hit_o.emtfTheta1() << std::endl << std::endl;
    std::cout << "WG  " << hit_o.wire1() << " - Phi: " << hit_o.emtfPhi() << std::endl;
  }
  std::cout << std::endl;
  }
 
  std::cout << std::endl;
}



