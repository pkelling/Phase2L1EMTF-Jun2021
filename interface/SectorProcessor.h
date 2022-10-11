#ifndef L1Trigger_Phase2L1EMTF_SectorProcessor_h
#define L1Trigger_Phase2L1EMTF_SectorProcessor_h

#include <type_traits>

#include "DataFormats/Provenance/interface/EventID.h"

#include "L1Trigger/Phase2L1EMTF/interface/Common.h"
#include "L1Trigger/Phase2L1EMTF/interface/SubsystemTags.h"
#include "L1Trigger/Phase2L1EMTF/interface/SubsystemCollection.h"

namespace emtf {

  namespace phase2 {

    class EMTFWorker;

    class SectorProcessor {
    public:
      void process(const EMTFWorker& iWorker,
                   int endcap,
                   int sector,
                   const edm::EventID& evt_id,
                   const SubsystemCollection& muon_primitives,
                   EMTFHitCollection& out_hits,
                   EMTFTrackCollection& out_tracks) const;

    private:
      template <typename>
      struct dependent_false;

      void process_step_1(const EMTFWorker& iWorker,
                          int endcap,
                          int sector,
                          int bx,
                          const SubsystemCollection& muon_primitives,
                          EMTFHitCollection& sector_hits) const;

      void process_step_2(const EMTFWorker& iWorker,
                          int endcap,
                          int sector,
                          int bx,
                          const EMTFHitCollection& sector_hits,
                          EMTFTrackCollection& sector_tracks) const;

      void dump_input_output(const edm::EventID& evt_id,
                             const SubsystemCollection& muon_primitives,
                             const EMTFHitCollection& out_hits,
                             const EMTFTrackCollection& out_tracks) const;


      int ge11_phi_conv(int chamber, int ph_init, int hs) const;
      int ge21_phi_conv(int chamber, int ph_init, int hs) const;
      void create_csc_luts(const EMTFWorker& iWorker) const;
      void create_gem_luts(const EMTFWorker& iWorker) const;
      void create_me0_luts(const EMTFWorker& iWorker) const;

      void test_dump(const EMTFWorker& iWorker, const SubsystemCollection& muon_primitives) const;

    };

    // Implementation of the templated classes and functions

    // A type-dependent expression that is always false
    template <typename>
    struct SectorProcessor::dependent_false : std::false_type {};

  }  // namespace phase2

}  // namespace emtf

#endif  // L1Trigger_Phase2L1EMTF_SectorProcessor_h not defined
