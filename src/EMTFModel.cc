#include "L1Trigger/Phase2L1EMTF/interface/EMTFModel.h"

// Xilinx HLS
#include "ap_int.h"
#include "ap_fixed.h"

// EMTF HLS
#include "emtf_hlslib.h"
#include "alternate.h"

using namespace emtf::phase2;

EMTFModel::EMTFModel(unsigned version, bool unconstrained) : version_(version), unconstrained_(unconstrained) {}

EMTFModel::~EMTFModel() {}

NdArrayDesc EMTFModel::get_input_shape() const {
  if (version_ == 3) {
    int n0 = num_emtf_chambers_v3 * num_emtf_segments_v3;
    int n1 = num_emtf_variables_v3;
    return NdArrayDesc({n0, n1});
  }
  return NdArrayDesc{};
}

NdArrayDesc EMTFModel::get_output_shape() const {
  if (version_ == 3) {
    int n0 = num_emtf_tracks_v3;
    int n1 = num_emtf_trk_variables_v3;
    return NdArrayDesc({n0, n1});
  }
  return NdArrayDesc{};
}

int EMTFModel::get_num_segments() const {
  if (version_ == 3) {
    return num_emtf_segments_v3;
  }
  return 0;
}

int EMTFModel::get_num_tracks() const {
  if (version_ == 3) {
    return num_emtf_tracks_v3;
  }
  return 0;
}

void EMTFModel::fit_impl_v3(const Vector& in0, Vector& out) const {
  // Check consistency with the parameters from namespace emtf_hlslib
  static_assert(EMTFModel::num_emtf_chambers_v3 == emtf_hlslib::phase2::num_emtf_chambers);
  static_assert(EMTFModel::num_emtf_segments_v3 == emtf_hlslib::phase2::num_emtf_segments);
  static_assert(EMTFModel::num_emtf_variables_v3 == emtf_hlslib::phase2::num_emtf_variables);
  static_assert(EMTFModel::num_emtf_tracks_v3 == emtf_hlslib::phase2::num_emtf_tracks);
  static_assert(EMTFModel::num_emtf_trk_variables_v3 ==
                (emtf_hlslib::phase2::num_emtf_features + emtf_hlslib::phase2::num_emtf_sites + 2));

  using namespace emtf_hlslib::phase2;

  // Unpack from in0
  // Note: the following are currently unused and will be synthesized away
  // - emtf_qual2, emtf_time, seg_cscfr, seg_gemdl, seg_bx
  emtf_phi_t emtf_phi[model_config::n_in];
  emtf_bend_t emtf_bend[model_config::n_in];
  emtf_theta1_t emtf_theta1[model_config::n_in];
  emtf_theta2_t emtf_theta2[model_config::n_in];
  emtf_qual1_t emtf_qual1[model_config::n_in];
  emtf_qual2_t emtf_qual2[model_config::n_in];
  emtf_time_t emtf_time[model_config::n_in];
  seg_zones_t seg_zones[model_config::n_in];
  seg_tzones_t seg_tzones[model_config::n_in];
  seg_cscfr_t seg_cscfr[model_config::n_in];
  seg_gemdl_t seg_gemdl[model_config::n_in];
  seg_bx_t seg_bx[model_config::n_in];
  seg_valid_t seg_valid[model_config::n_in];

  // Loop over in0
  auto in0_iter = in0.begin();

  for (unsigned iseg = 0; iseg < model_config::n_in; iseg++) {
    emtf_phi[iseg] = *(in0_iter++);
    emtf_bend[iseg] = *(in0_iter++);
    emtf_theta1[iseg] = *(in0_iter++);
    emtf_theta2[iseg] = *(in0_iter++);
    emtf_qual1[iseg] = *(in0_iter++);
    emtf_qual2[iseg] = *(in0_iter++);
    emtf_time[iseg] = *(in0_iter++);
    seg_zones[iseg] = *(in0_iter++);
    seg_tzones[iseg] = *(in0_iter++);
    seg_cscfr[iseg] = *(in0_iter++);
    seg_gemdl[iseg] = *(in0_iter++);
    seg_bx[iseg] = *(in0_iter++);
    seg_valid[iseg] = *(in0_iter++);
  }  // end loop over in0

  // Intermediate arrays (for layers 0..3)
  zoning_out_t zoning_0_out[zoning_config::n_out];
  zoning_out_t zoning_1_out[zoning_config::n_out];
  zoning_out_t zoning_2_out[zoning_config::n_out];
  pooling_out_t pooling_0_out[pooling_config::n_out];
  pooling_out_t pooling_1_out[pooling_config::n_out];
  pooling_out_t pooling_2_out[pooling_config::n_out];
  zonesorting_out_t zonesorting_0_out[zonesorting_config::n_out];
  zonesorting_out_t zonesorting_1_out[zonesorting_config::n_out];
  zonesorting_out_t zonesorting_2_out[zonesorting_config::n_out];
  zonemerging_out_t zonemerging_0_out[zonemerging_config::n_out];

  // Layer 0 - Zoning

  zoning_layer<m_zone_any_tag>(emtf_phi, seg_zones, seg_tzones, seg_valid, zoning_0_out, zoning_1_out, zoning_2_out);

  // Layer 1 - Pooling

  pooling_layer<m_zone_0_tag>(zoning_0_out, pooling_0_out);
  pooling_layer<m_zone_1_tag>(zoning_1_out, pooling_1_out);
  pooling_layer<m_zone_2_tag>(zoning_2_out, pooling_2_out);

  // Layer 2 - Zone sorting

  zonesorting_layer<m_zone_any_tag>(pooling_0_out, zonesorting_0_out);
  zonesorting_layer<m_zone_any_tag>(pooling_1_out, zonesorting_1_out);
  zonesorting_layer<m_zone_any_tag>(pooling_2_out, zonesorting_2_out);

  // Layer 3 - Zone merging

  zonemerging_layer<m_zone_any_tag>(zonesorting_0_out, zonesorting_1_out, zonesorting_2_out, zonemerging_0_out);

  // Unpack from in1 (a.k.a. zonemerging_0_out)
  trk_qual_t trk_qual[trkbuilding_config::n_in];
  trk_patt_t trk_patt[trkbuilding_config::n_in];
  trk_col_t trk_col[trkbuilding_config::n_in];
  trk_zone_t trk_zone[trkbuilding_config::n_in];
  trk_tzone_t trk_tzone[trkbuilding_config::n_in];

  // Loop over in1
  for (unsigned itrk = 0; itrk < trkbuilding_config::n_in; itrk++) {
    const trkbuilding_in_t curr_trk_in = zonemerging_0_out[itrk];
    const trk_tzone_t curr_trk_tzone = detail::timezone_traits<m_timezone_0_tag>::value;  // default timezone

    constexpr int bits_lo_0 = 0;
    constexpr int bits_lo_1 = trk_qual_t::width;
    constexpr int bits_lo_2 = pooling_out_t::width;
    constexpr int bits_lo_3 = zonesorting_out_t::width;
    constexpr int bits_lo_4 = zonemerging_out_t::width;

    trk_qual[itrk] = curr_trk_in.range(bits_lo_1 - 1, bits_lo_0);
    trk_patt[itrk] = curr_trk_in.range(bits_lo_2 - 1, bits_lo_1);
    trk_col[itrk] = curr_trk_in.range(bits_lo_3 - 1, bits_lo_2);
    trk_zone[itrk] = curr_trk_in.range(bits_lo_4 - 1, bits_lo_3);
    trk_tzone[itrk] = curr_trk_tzone;
  }  // end loop over in1

  // Intermediate arrays (for layers 4..6)
  trk_seg_t trk_seg[trkbuilding_config::n_out * num_emtf_sites];
  trk_seg_v_t trk_seg_v[trkbuilding_config::n_out];
  trk_feat_t trk_feat[trkbuilding_config::n_out * num_emtf_features];
  trk_valid_t trk_valid[trkbuilding_config::n_out];
  trk_seg_t trk_seg_rm[duperemoval_config::n_out * num_emtf_sites];
  trk_seg_v_t trk_seg_rm_v[duperemoval_config::n_out];
  trk_feat_t trk_feat_rm[duperemoval_config::n_out * num_emtf_features];
  trk_valid_t trk_valid_rm[duperemoval_config::n_out];
  trk_origin_t trk_origin_rm[duperemoval_config::n_out];
  trk_invpt_t trk_invpt[fullyconnect_config::n_out];
  trk_phi_t trk_phi[fullyconnect_config::n_out];
  trk_eta_t trk_eta[fullyconnect_config::n_out];
  trk_d0_t trk_d0[fullyconnect_config::n_out];
  trk_z0_t trk_z0[fullyconnect_config::n_out];
  trk_beta_t trk_beta[fullyconnect_config::n_out];

  // Layer 4 - Track building

  for (unsigned itrk = 0; itrk < trkbuilding_config::n_in; itrk++) {
    // Intermediate arrays (for layer output)
    trk_seg_t curr_trk_seg[num_emtf_sites];
    trk_feat_t curr_trk_feat[num_emtf_features];

    trkbuilding_layer<m_zone_any_tag>(emtf_phi,
                                      emtf_bend,
                                      emtf_theta1,
                                      emtf_theta2,
                                      emtf_qual1,
                                      emtf_qual2,
                                      emtf_time,
                                      seg_zones,
                                      seg_tzones,
                                      seg_cscfr,
                                      seg_gemdl,
                                      seg_bx,
                                      seg_valid,
                                      trk_qual[itrk],
                                      trk_patt[itrk],
                                      trk_col[itrk],
                                      trk_zone[itrk],
                                      trk_tzone[itrk],
                                      curr_trk_seg,
                                      trk_seg_v[itrk],
                                      curr_trk_feat,
                                      trk_valid[itrk]);

    // Copy to arrays
    detail::copy_n_values<num_emtf_sites>(curr_trk_seg, &(trk_seg[itrk * num_emtf_sites]));
    detail::copy_n_values<num_emtf_features>(curr_trk_feat, &(trk_feat[itrk * num_emtf_features]));
  }  // end loop over tracks

  // Layer 5 - Duplicate removal

  duperemoval_layer<m_zone_any_tag>(
      trk_seg, trk_seg_v, trk_feat, trk_valid, trk_seg_rm, trk_seg_rm_v, trk_feat_rm, trk_valid_rm, trk_origin_rm);

  // Layer 6 - Fully connected

  auto fullyconnect_layer_cheat = [](const trk_feat_t curr_trk_feat_rm[num_emtf_features],
                                     trk_invpt_t& curr_trk_invpt,
                                     trk_phi_t& curr_trk_phi,
                                     trk_eta_t& curr_trk_eta,
                                     trk_d0_t& curr_trk_d0,
                                     trk_z0_t& curr_trk_z0,
                                     trk_beta_t& curr_trk_beta) -> void {
    curr_trk_invpt.range() = trk_invpt_default_v3;
    curr_trk_phi.range() = trk_phi_default_v3;
    curr_trk_eta.range() = trk_eta_default_v3;
    curr_trk_d0.range() = trk_d0_default_v3;
    curr_trk_z0.range() = trk_z0_default_v3;
    curr_trk_beta.range() = trk_beta_default_v3;
  };

  for (unsigned itrk = 0; itrk < fullyconnect_config::n_in; itrk++) {
    // Skip fullyconnect_layer if invalid
    if (not trk_valid_rm[itrk]) {
      fullyconnect_layer_cheat(
          nullptr, trk_invpt[itrk], trk_phi[itrk], trk_eta[itrk], trk_d0[itrk], trk_z0[itrk], trk_beta[itrk]);
      continue;
    }

    // Intermediate arrays (for layer input)
    trk_feat_t curr_trk_feat_rm[num_emtf_features];

    // Copy from arrays
    detail::copy_n_values<num_emtf_features>(&(trk_feat_rm[itrk * num_emtf_features]), curr_trk_feat_rm);

    fullyconnect_layer<m_zone_any_tag>(
        curr_trk_feat_rm, trk_invpt[itrk], trk_phi[itrk], trk_eta[itrk], trk_d0[itrk], trk_z0[itrk], trk_beta[itrk]);
  }  // end loop over tracks

  // Copy to output: trk_feat_rm, trk_seg_rm, trk_valid_rm, trk_invpt
  auto out_iter = out.begin();

  for (unsigned i = 0; i < model_config::n_out; i++) {
    const unsigned itrk = (i / model_config::n_out_per_trk);
    const unsigned ivar = (i % model_config::n_out_per_trk);

    auto curr_trk_seg_rm = &(trk_seg_rm[itrk * num_emtf_sites]);
    auto curr_trk_feat_rm = &(trk_feat_rm[itrk * num_emtf_features]);

    if (ivar < num_emtf_features) {
      *(out_iter++) = curr_trk_feat_rm[ivar];
    } else if (ivar < (num_emtf_features + num_emtf_sites + 0)) {
      const unsigned ivar_1 = (ivar - num_emtf_features);
      const trk_seg_t invalid_marker_trk_seg = model_config::n_in;
      *(out_iter++) = (trk_seg_rm_v[itrk][ivar_1]) ? curr_trk_seg_rm[ivar_1] : invalid_marker_trk_seg;
    } else if (ivar < (num_emtf_features + num_emtf_sites + 1)) {
      *(out_iter++) = trk_valid_rm[itrk];
    } else if (ivar < (num_emtf_features + num_emtf_sites + 2)) {
      const trk_invpt_t invalid_marker_trk_invpt = ap_int_limits<trk_invpt_t>::min_value;
      *(out_iter++) = (trk_valid_rm[itrk]) ? trk_invpt[itrk] : invalid_marker_trk_invpt;
    }
  }  // end loop over out


      /************** Alternate inputs by site (w/new definition of site - ME & RE ring 1 and 2 are separate) *****************/

      prompt_trkbuild::emtf_seg_t emtf_segs[prompt_trkbuild::num_sites][prompt_trkbuild::max_ch_site][prompt_trkbuild::max_seg_ch];

      prompt_trkbuild::convert_input(emtf_phi, emtf_theta1, emtf_theta2, emtf_bend,
                                     emtf_qual1, emtf_time, seg_zones, seg_tzones,
                                     seg_bx, seg_valid, emtf_segs);

      prompt_trkbuild::hitmap_t hitmaps[3][8];
      prompt_trkbuild::build_hitmaps(emtf_segs, hitmaps);

      prompt_trkbuild::patt_match_t patt_matches[3][288];
      prompt_trkbuild::match_patterns(hitmaps, patt_matches);

      prompt_trkbuild::best_trk_t best_tracks[4];
      prompt_trkbuild::find_best_tracks(patt_matches, best_tracks);

      std::cout << std::endl;
      prompt_trkbuild::trk_feat_t trk_features_alt[4][40]; 
      prompt_trkbuild::seg_id_t trk_seg_ids_alt[4][12];
      prompt_trkbuild::trk_building(best_tracks, emtf_segs, trk_features_alt, trk_seg_ids_alt);
      prompt_trkbuild::trk_building(best_tracks, emtf_segs, trk_features_test, trk_seg_ids, ph_median_test, th_median_test);
      /*********************************************************************************************************************/


      /****************************************** TEST ZONING *****************************************************************/
      int passed = 1;
      for(int zone=0; zone<3; zone++){
          for(int row=0; row<8; row++){

            ap_uint<288> target;
              if( zone == 0)
                  target = zoning_0_out[row];
              else if(zone == 1)
                  target = zoning_1_out[row];
              else
                  target = zoning_2_out[row];

              if(target != hitmaps[zone][row]){
                  std::cout << "\n\nFailed Comparison. \tZone: " << zone << "\tRow: " << row << std::endl;

                  std::cout << std::hex << target << std::endl;
                  std::cout << std::hex << hitmaps[zone][row] << "\n\n" << std::dec;
                  passed = 0;
              }
          }
      }

      /****************************************** TEST POOLING *****************************************************************/
      for(int zone=0; zone<3; zone++){
          for(int col=0; col<288; col++){
            ap_uint<9> target;
              if( zone == 0)
                  target = pooling_0_out[col];
              else if(zone == 1)
                  target = pooling_1_out[col];
              else
                  target = pooling_2_out[col];

              if(target != ap_uint<9>( (patt_matches[zone][col].patt,patt_matches[zone][col].qual)) ){
                  std::cout << "\n\nFailed Comparison. \tZone: " << zone << "\tColumn: " << col << std::endl;
                  std::cout << std::hex << target << std::endl;
                  std::cout << std::hex << ap_uint<9>((patt_matches[zone][col].patt,patt_matches[zone][col].qual)) << "\n\n" << std::dec;
                  passed = 0;
              }
          }
      }

      /****************************************** TEST SORTING AND MERGING *****************************************************************/
      for(int trk=0; trk<4; trk++){
        if( zonemerging_0_out[trk].range(5,0) != 0 and zonemerging_0_out[trk] != ap_uint<20>((best_tracks[trk].zone, best_tracks[trk].col, best_tracks[trk].patt, best_tracks[trk].qual)) ){
            std::cout << "\n\nFailed Comparison. \tTrack: " << trk << std::endl;
            std::cout << "Quality: " << zonemerging_0_out[trk].range(5,0) << " vs " << best_tracks[trk].qual << std::endl;
            std::cout << "Pattern: " << zonemerging_0_out[trk].range(8,6) << " vs " << best_tracks[trk].patt << std::endl;
            std::cout << "Column: " << zonemerging_0_out[trk].range(17,9) << " vs " << best_tracks[trk].col << std::endl;
            std::cout << "Zone: " << zonemerging_0_out[trk].range(19,18) << " vs " << best_tracks[trk].zone << std::endl;
        }
      }

      /****************************************** TEST Track Building ***********************************************************************/

      for(int trk=0; trk<4; trk++){
          for(int feat=0; feat<40; feat++){
              if( trk_features_alt[trk][feat] != trk_feat[trk *num_emtf_features + feat]){
                  std::cout << trk_features_alt[trk][feat] << " vs " << trk_feat[trk *num_emtf_features + feat] << "\t";
                  passed = 0;
              }
          }
      }



      if(!passed){
          for (unsigned iseg = 0; iseg < model_config::n_in; iseg++) {
              if(seg_valid[iseg]){
                  std::cout << "\n" << "Seg: " << iseg << " - ";
                  std::cout << emtf_phi[iseg] << ", ";
                  std::cout << emtf_bend[iseg]  << ", ";
                  std::cout << emtf_theta1[iseg]  << ", ";
                  std::cout << emtf_theta2[iseg]  << ", ";
                  std::cout << emtf_qual1[iseg]  << ", ";
                  std::cout << emtf_qual2[iseg]  << ", ";
                  std::cout << emtf_time[iseg]  << ", ";
                  std::cout << seg_zones[iseg]  << ", ";
                  std::cout << seg_tzones[iseg]  << ", ";
                  std::cout << seg_cscfr[iseg]  << ", ";
                  std::cout << seg_gemdl[iseg]  << ", ";
                  std::cout << seg_bx[iseg]  << ", ";
                  std::cout << seg_valid[iseg];
              }
          }  // end loop over in0
          std::cout << std::endl;
      }

  /*********************************************************************************************************************/

}
