#include "L1Trigger/Phase2L1EMTF/test/tools/Toolbox.h"

namespace emtf {

  namespace phase2 {

    namespace toolbox {

      int get_csc_chamber(int station, int sector, int subsector, int ring, int cscid) {
        constexpr int kInvalid = -99;
        int chamber = kInvalid;
        if (station == 1) {
          // Chamber offset of 2: First chamber in sector 1 is chamber 3
          chamber = ((sector - 1) * 6) + cscid + 2;
          if (ring == 2)
            chamber -= 3;
          else if (ring == 3)
            chamber -= 6;
          if (subsector == 2)
            chamber += 3;
          if (chamber > 36)
            chamber -= 36;
        } else {
          if (ring == 1) {
            // Chamber offset of 1: First chamber in sector 1 is chamber 2
            chamber = ((sector - 1) * 3) + cscid + 1;
            if (chamber > 18)
              chamber -= 18;
          } else {
            // Chamber offset of 2: First chamber in sector 1 is chamber 3
            chamber = ((sector - 1) * 6) + cscid - 3 + 2;
            if (chamber > 36)
              chamber -= 36;
          }
        }
        return chamber;
      }

      int next_csc_chamber_10deg(int chamber) { return (chamber == 36) ? 1 : chamber + 1; }

      int prev_csc_chamber_10deg(int chamber) { return (chamber == 1) ? 36 : chamber - 1; }

      int next_csc_chamber_20deg(int chamber) { return (chamber == 18) ? 1 : chamber + 1; }

      int prev_csc_chamber_20deg(int chamber) { return (chamber == 1) ? 18 : chamber - 1; }

      int get_csc_ring(int station, int cscid, int strip) {
        constexpr int kInvalid = -99;
        if (station == 1) {
          if (cscid < 4 && strip >= 128)
            return 4;
          else if (cscid < 4)
            return 1;
          else if (cscid < 7)
            return 2;
          else if (cscid < 10)
            return 3;
          else
            return kInvalid;
        } else {
          if (cscid < 4)
            return 1;
          else if (cscid < 10)
            return 2;
          else
            return kInvalid;
        }
      }

      int get_trigger_subsector(int station, int chamber) { return (station != 1) ? 0 : (((chamber % 6) > 2) ? 1 : 2); }

      // Copied from DataFormats/MuonDetId/src/CSCDetId.cc
      int get_trigger_sector(int ring, int station, int chamber) {
        int result = 0;
        if (station > 1 && ring > 1) {
          result = ((static_cast<unsigned>(chamber - 3) & 0x7f) / 6) + 1;  // ch 3-8->1, 9-14->2, ... 1,2 -> 6
        } else if (station == 1) {
          result = ((static_cast<unsigned>(chamber - 3) & 0x7f) / 6) + 1;  // ch 3-8->1, 9-14->2, ... 1,2 -> 6
        } else {
          result = ((static_cast<unsigned>(chamber - 2) & 0x1f) / 3) + 1;  // ch 2-4-> 1, 5-7->2, ...
        }
        // max sector is 6, some calculations give a value greater than six but this is expected.
        return (result <= 6) ? result : 6;
      }

      int next_trigger_sector(int sector) { return (sector == 6) ? 1 : sector + 1; }

      int prev_trigger_sector(int sector) { return (sector == 1) ? 6 : sector - 1; }

      // Copied from DataFormats/MuonDetId/src/CSCDetId.cc
      int get_trigger_cscid(int ring, int station, int chamber) {
        int result = 0;
        if (station == 1) {
          result = (chamber) % 3 + 1;  // 1,2,3
          switch (ring) {
            case 1:
            case 4:
              break;
            case 2:
              result += 3;  // 4,5,6
              break;
            case 3:
              result += 6;  // 7,8,9
              break;
          }
        } else {
          if (ring == 1) {
            result = (chamber + 1) % 3 + 1;  // 1,2,3
          } else {
            result = (chamber + 3) % 6 + 4;  // 4,5,6,7,8,9
          }
        }
        return result;
      }

      // Copied from RecoMuon/DetLayers/src/MuonCSCDetLayerGeometryBuilder.cc
      bool get_trigger_cscfr(int ring, int station, int chamber) {
        bool result = false;

        bool isOverlapping = !(station == 1 && ring == 3);
        // not overlapping means back
        if (isOverlapping) {
          bool isEven = (chamber % 2 == 0);
          // odd chambers are bolted to the iron, which faces
          // forward in 1&2, backward in 3&4, so...
          result = (station < 3) ? isEven : !isEven;
        }
        return result;
      }

      // Numbers of halfstrips and wiregroups
      //
      // +----------------------------+------------+------------+
      // | Chamber type               | Num of     | Num of     |
      // |                            | halfstrips | wiregroups |
      // +----------------------------+------------+------------+
      // | ME1/1a                     | 96         | 48         |
      // | ME1/1b                     | 128        | 48         |
      // | ME1/2                      | 160        | 64         |
      // | ME1/3                      | 128        | 32         |
      // | ME2/1                      | 160        | 112        |
      // | ME3/1, ME4/1               | 160        | 96         |
      // | ME2/2, ME3/2, ME4/2        | 160        | 64         |
      // +----------------------------+------------+------------+
      std::pair<int, int> get_csc_max_strip_and_wire(int station, int ring) {
        int max_strip = 0;                // halfstrip
        int max_wire = 0;                 // wiregroup
        if (station == 1 && ring == 4) {  // ME1/1a
          max_strip = 96;
          max_wire = 48;
        } else if (station == 1 && ring == 1) {  // ME1/1b
          max_strip = 128;
          max_wire = 48;
        } else if (station == 1 && ring == 2) {  // ME1/2
          max_strip = 160;
          max_wire = 64;
        } else if (station == 1 && ring == 3) {  // ME1/3
          max_strip = 128;
          max_wire = 32;
        } else if (station == 2 && ring == 1) {  // ME2/1
          max_strip = 160;
          max_wire = 112;
        } else if (station >= 3 && ring == 1) {  // ME3/1, ME4/1
          max_strip = 160;
          max_wire = 96;
        } else if (station >= 2 && ring == 2) {  // ME2/2, ME3/2, ME4/2
          max_strip = 160;
          max_wire = 64;
        }
        return std::make_pair(max_strip, max_wire);
      }

      // Numbers of patterns and qualities
      std::pair<int, int> get_csc_max_pattern_and_quality(int station, int ring) {
        int max_pattern = 11;
        int max_quality = 16;
        return std::make_pair(max_pattern, max_quality);
      }

      // Calculates special chamber ID for track address sent to uGMT by using cscid, subsector, neighbor, and station
      // Please refers to DN-2015/017 for uGMT conventions
      int get_ugmt_chamber(int cscid, int subsector, int neighbor, int station) {
        if (station == 1) {
          if (cscid == 3 && neighbor == 1 && subsector == 2)
            return 1;
          else if (cscid == 6 && neighbor == 1 && subsector == 2)
            return 2;
          else if (cscid == 9 && neighbor == 1 && subsector == 2)
            return 3;
          else if (cscid == 3 && neighbor == 0 && subsector == 2)
            return 4;
          else if (cscid == 6 && neighbor == 0 && subsector == 2)
            return 5;
          else if (cscid == 9 && neighbor == 0 && subsector == 2)
            return 6;
          else
            return 0;
        } else {
          if (cscid == 3 && neighbor == 1)
            return 1;
          else if (cscid == 9 && neighbor == 1)
            return 2;
          else if (cscid == 3 && neighbor == 0)
            return 3;
          else if (cscid == 9 && neighbor == 0)
            return 4;
          else
            return 0;
        }
      }

    }  // namespace toolbox

  }  // namespace phase2

}  // namespace emtf
