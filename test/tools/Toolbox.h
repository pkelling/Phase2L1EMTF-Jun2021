#ifndef L1Trigger_Phase2L1EMTF_Toolbox_h
#define L1Trigger_Phase2L1EMTF_Toolbox_h

#include <cmath>
#include <utility>  // provides std::pair

namespace emtf {

  namespace phase2 {

    namespace toolbox {

      // CSCDetId chamber number
      int get_csc_chamber(int station, int sector, int subsector, int ring, int cscid);

      int next_csc_chamber_10deg(int chamber);

      int prev_csc_chamber_10deg(int chamber);

      int next_csc_chamber_20deg(int chamber);

      int prev_csc_chamber_20deg(int chamber);

      // CSCDetId ring number (assume ME1/1a -> ring 4 convention)
      int get_csc_ring(int station, int cscid, int strip);

      // CSC trigger subsector
      int get_trigger_subsector(int station, int chamber);

      // CSC trigger sector
      int get_trigger_sector(int ring, int station, int chamber);

      int next_trigger_sector(int sector);

      int prev_trigger_sector(int sector);

      // CSC trigger cscid
      int get_trigger_cscid(int ring, int station, int chamber);

      // CSC trigger front/rear bit
      bool get_trigger_cscfr(int ring, int station, int chamber);

      // CSC max strip & max wire
      std::pair<int, int> get_csc_max_strip_and_wire(int station, int ring);

      // CSC max pattern & max quality
      std::pair<int, int> get_csc_max_pattern_and_quality(int station, int ring);

      // uGMT chamber number
      int get_ugmt_chamber(int cscid, int subsector, int neighbor, int station);

      // _______________________________________________________________________
      // radians <-> degrees
      inline constexpr float deg_to_rad(float deg) {
        constexpr float factor = M_PI / 180.;
        return deg * factor;
      }

      inline constexpr float rad_to_deg(float rad) {
        constexpr float factor = 180. / M_PI;
        return rad * factor;
      }

      // _______________________________________________________________________
      // phi range: [-180..180] or [-pi..pi]
      inline constexpr float wrap_phi_deg(float deg) {
        constexpr float twopi = 360.;
        constexpr float recip = 1.0 / twopi;
        return deg - (std::round(deg * recip) * twopi);
      }

      inline constexpr float wrap_phi_rad(float rad) {
        constexpr float twopi = M_PI * 2.;
        constexpr float recip = 1.0 / twopi;
        return rad - (std::round(rad * recip) * twopi);
      }

      // _______________________________________________________________________
      // eta
      inline float calc_eta_from_theta_rad(float theta_rad) {
        float eta = -1. * std::log(std::tan(theta_rad / 2.));
        return eta;
      }

      inline float calc_eta_from_theta_deg(float theta_deg) {
        float eta = calc_eta_from_theta_rad(deg_to_rad(theta_deg));
        return eta;
      }

      // _______________________________________________________________________
      // theta
      inline float calc_theta_rad_from_eta(float eta) {
        float theta = std::atan2(1.0, std::sinh(eta));  // cot(theta) = sinh(eta)
        return theta;
      }

      inline float calc_theta_deg_from_eta(float eta) {
        float theta = rad_to_deg(calc_theta_rad_from_eta(eta));
        return theta;
      }

      inline float calc_theta_deg_from_int(int theta_int) {
        float theta = static_cast<float>(theta_int);
        theta = theta * (45.0 - 8.5) / 128. + 8.5;
        return theta;
      }

      inline float calc_theta_rad_from_int(int theta_int) {
        float theta = deg_to_rad(calc_theta_deg_from_int(theta_int));
        return theta;
      }

      inline int calc_theta_int(float theta, int endcap) {  // theta in deg [0..180], endcap [-1, +1]
        theta = (endcap == -1) ? (180. - theta) : theta;
        theta = (theta - 8.5) * 128. / (45.0 - 8.5);
        int theta_int = static_cast<int>(std::round(theta));
        theta_int = (theta_int <= 0) ? 1 : theta_int;  // protect against invalid value
        return theta_int;
      }

      // _______________________________________________________________________
      // phi
      inline float calc_phi_glob_deg_from_loc(float loc, int sector) {  // loc in deg, sector [1..6]
        float glob = loc + 15. + (60. * (sector - 1));
        glob = (glob >= 180.) ? (glob - 360.) : glob;
        return glob;
      }

      inline float calc_phi_glob_rad_from_loc(float loc, int sector) {  // loc in rad, sector [1..6]
        float glob = deg_to_rad(calc_phi_glob_deg_from_loc(rad_to_deg(loc), sector));
        return glob;
      }

      inline float calc_phi_loc_deg_from_int(int phi_int) {
        float loc = static_cast<float>(phi_int);
        loc = (loc / 60.) - 22.;
        return loc;
      }

      inline float calc_phi_loc_rad_from_int(int phi_int) {
        float loc = deg_to_rad(calc_phi_loc_deg_from_int(phi_int));
        return loc;
      }

      inline float calc_phi_loc_deg_from_glob(float glob, int sector) {  // glob in deg [-180..180], sector [1..6]
        glob = wrap_phi_deg(glob);
        float loc = glob - 15. - (60. * (sector - 1));
        return loc;
      }

      inline int calc_phi_int(float glob, int sector) {  // glob in deg [-180..180], sector [1..6]
        float loc = calc_phi_loc_deg_from_glob(glob, sector);
        loc = ((loc + 22.) < 0.) ? (loc + 360.) : loc;
        loc = (loc + 22.) * 60.;
        int phi_int = static_cast<int>(std::round(loc));
        return phi_int;
      }

    }  // namespace toolbox

  }  // namespace phase2

}  // namespace emtf

#endif  // L1Trigger_Phase2L1EMTF_Toolbox_h not defined
