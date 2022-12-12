#include <cmath>
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"

#include "L1Trigger/Phase2L1EMTF/interface/Toolbox.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

using namespace emtf::phase2;


class MakePh2CoordLUT : public edm::EDAnalyzer {
public:
  explicit MakePh2CoordLUT(const edm::ParameterSet&);
  virtual ~MakePh2CoordLUT();

private:
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);

  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // Generate LUTs
  void generateLUTs();

  // CSC LUTs
  void generate_csc_LUTs(int endcap, int sector);
  void generate_ME11_csc_LUT(int endcap, int sector, int ch, int ich);
  void generate_nonME11_csc_LUT(int endcap, int sector, int st, int ring, int ch, int ich);
  GlobalPoint get_csc_global_point(const CSCLayerGeometry* layer_geom, const CSCLayer* layer, const int wiregroup, const int eighth_strip);

  // GEM/ME0 LUTs
  void generate_gem_LUTs(int endcap, int sector);
  void generate_me0_LUTs(int endcap, int sector);

  // Quick LUT Verification (just checks against hardcoded ph_init values)
  void quick_LUT_verification(int endcap, int sector);

private:
  const edm::ParameterSet config_;

  bool done_;

  /// Event setup
  const CSCGeometry* theCSCGeometry_;
  const GEMGeometry* theGEMGeometry_;
  const ME0Geometry* theME0Geometry_;
};

// _____________________________________________________________________________

MakePh2CoordLUT::MakePh2CoordLUT(const edm::ParameterSet& iConfig)
    : config_(iConfig),
      done_(false){
  std::cout << "Initializing MakePh2CoordLUT\n";
}

MakePh2CoordLUT::~MakePh2CoordLUT() {}


void MakePh2CoordLUT::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  std::cout << "Beginning MakePh2CoordLUT Run\n";

  /// Setup Geometries

  // CSC Geometry
  edm::ESHandle<CSCGeometry> cscGeometryHandle;
  iSetup.get<MuonGeometryRecord>().get(cscGeometryHandle);
  assert(cscGeometryHandle.isValid());
  theCSCGeometry_ = cscGeometryHandle.product();

  // GEM Geometry
  edm::ESHandle<GEMGeometry> gemGeometryHandle;
  iSetup.get<MuonGeometryRecord>().get(gemGeometryHandle);
  assert(gemGeometryHandle.isValid());
  theGEMGeometry_ = gemGeometryHandle.product();

  // ME0 Geometry
  edm::ESHandle<ME0Geometry> me0GeometryHandle;
  iSetup.get<MuonGeometryRecord>().get(me0GeometryHandle);
  assert(me0GeometryHandle.isValid());
  theME0Geometry_ = me0GeometryHandle.product();
}


void MakePh2CoordLUT::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  std::cout << "Ending MakePh2CoordLUT Run\n";
}


void MakePh2CoordLUT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (done_)
    return;
  
  std::cout << "Running Analyzer for MakePh2CoordLUT\n";
  generateLUTs();

  done_ = true;
  return;
}



// _____________ Generate the Lookup Tables -------------------
void MakePh2CoordLUT::generateLUTs() {
  std::cout << "Running generateLUTs()\n";

  int endcap = 1; // [+1,-1]
  int sector = 1; // [1..6]
  
  generate_csc_LUTs(endcap, sector); // endcap [+1,-1], sector [1..6]
  generate_me0_LUTs(endcap, sector); // endcap [+1,-1], sector [1..6]
  generate_gem_LUTs(endcap, sector); // endcap [+1,-1], sector [1..6]

  quick_LUT_verification(endcap, sector);
}


void MakePh2CoordLUT::generate_gem_LUTs(int endcap, int sector){
  // Notes:
  // - Creates LUT looping through LUT addresses and building the trigger primitives backwards from the address (bits dependent on station)
  //
  // Phi Values:
  // - Currently uses bottom of eta range for phi-address (GEM Geometry doesn't provide function to find sub-partition points)
  // - Halfpads either mean the center of a pad (even hps) or the middle of 2 pads (odd hps) - Note: GEMs have strips, but EMTF recieves pad info only

  int sector_idx = sector-1; // For start and end chs
  int first_sector_chambers_ge11[6] = {2, 8, 14, 20, 26, 32};
  int first_sector_chambers_ge21[6] = {1, 4, 7,  10, 13, 16};

  for(int st=1; st<=2; st++){ // ge11 and ge21
    // Station Independent variables
    int num_theta_bits = 3; // only 8 partitions (numbered 1->8 in cmssw)

    // Station Dependent Variables
    int num_chambers  = (st == 1)? 7 : 4; // ge11 is 10degs, ge21 is 20degs
    int chamber       = (st == 1)? first_sector_chambers_ge11[sector_idx] : first_sector_chambers_ge21[sector_idx];
    int num_halfpads  = (st == 1)? (192*2-1)  : (384*2-1); // pads: ge11-192, ge21-384 :  -1 because last halfpad is the center of the last pad (not it's edge)
    int num_hp_bits   = (st == 1)? 8 : 9; // could change this to 8 for ge11
    int phi_addr_eta_bits = (st == 1)? 3 : 2;

    for(int ich=0; ich<num_chambers; ich++){ 
      // Get real chamber number
      chamber += ich;
      if(st==1 && chamber > 36) 
        chamber -= 36;
      else if(st==2 && chamber > 18)
        chamber -= 18;

      // GEM Phi LUT
      int num_phi_addresses = std::pow(2,num_hp_bits + phi_addr_eta_bits);
      std::vector<int> chamber_phi_LUT; 
      chamber_phi_LUT.reserve(num_phi_addresses); // Need to fill whole address space

      for(int addr=0; addr<num_phi_addresses; addr++){
        int phi_int = 0;

        // Get TP info from address
        int halfpad = (addr >> phi_addr_eta_bits);
        int roll_bitmask = (1 << phi_addr_eta_bits) - 1; // masks bits used as eta info (0x7 or 0x3)
        int iroll = ((addr & roll_bitmask) << (num_theta_bits - phi_addr_eta_bits)) + 1; // Gets tp partition from address (lowest of possible ones) [1 indexed]

        if(halfpad < num_halfpads){
          GEMDetId detid(endcap, 1, st, 1, chamber, iroll); // (region, ring, station, layer, chamber, roll )
          const GEMEtaPartition* roll = theGEMGeometry_->etaPartition(detid);

          const float center_of_pad = 0.5f + (0.5f * halfpad); // halfpad=0 is the center of the first pad
          const LocalPoint& lp = roll->centreOfPad(center_of_pad);
          const GlobalPoint& gp = roll->surface().toGlobal(lp);

          float phi_degs = toolbox::rad_to_deg(gp.barePhi());
          phi_int = toolbox::calc_phi_int(phi_degs, sector);
        }

        chamber_phi_LUT.push_back(phi_int);
      }


      // GEM Theta LUT
      int num_theta_addresses = std::pow(2,num_theta_bits);
      std::vector<int> chamber_theta_LUT; 
      chamber_theta_LUT.reserve(num_theta_addresses); // Need to fill whole addres space

      for(int addr=0; addr<num_theta_addresses; addr++){
        int iroll = addr + 1; // cmssw rolls are 1 indexed
        int halfpad = num_halfpads / 2; // use middle of chamber for eta

        GEMDetId detid(endcap, 1, st, 1, chamber, iroll); // (region, ring, station, layer, chamber, roll )
        const GEMEtaPartition* roll = theGEMGeometry_->etaPartition(detid);

        const float center_of_pad = 0.5f + (0.5f * halfpad); 
        const LocalPoint& lp = roll->centreOfPad(center_of_pad);
        const GlobalPoint& gp = roll->surface().toGlobal(lp);

        float theta_degs = toolbox::rad_to_deg(gp.theta());
        int theta_int = toolbox::calc_theta_int(theta_degs, endcap);
        chamber_theta_LUT.push_back(theta_int);
      }

      // Fill GEM Phi LUT
      std::stringstream phi_filename; 
      phi_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/GE" << st << "1_phi_ch" << ich << ".txt";
      std::ofstream phiFile(phi_filename.str());
      for (auto phi : chamber_phi_LUT) phiFile << phi << ", ";
      phiFile.close();

      // Fill GEM Theta LUT
      std::stringstream theta_filename; 
      theta_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/GE" << st << "1_theta_ch" << ich << ".txt";
      std::ofstream thetaFile(theta_filename.str());
      for (auto theta : chamber_theta_LUT) thetaFile << theta << ", ";
      thetaFile.close();

    }
  }

}



void MakePh2CoordLUT::generate_me0_LUTs(int endcap, int sector){
  // Notes:
  // - Creates LUT looping through LUT addresses and building the trigger primitives backwards from the address {10bit hs, 3bit eta}
  // -eta number (4 bits) represent either a single partition hit (even number) or a hit that spans 2 partitions (odd number)
  //      - Currenly, ME0 Geometry can only get theta for a partition, you cannot get a sub-partition
  //      - I dont think Eta number 15 should ever occur
  //      - If you decrease eta_bits in phi address, you should change tp_partition to middle of range
  //      - Also, I think using centreOfStrip() gives the phi position based on the bottom of the strip (bottom of partition, vs the middle)

  int num_halfstrips = 768; // 2 * n_me0_strips
  int num_hs_bits = 10;
  int num_eta_bits = 4; // gives half-partition resolution
  int phi_addr_eta_bits = 3;

  int first_sector_chamber[6] = {1,4,7,10,13,16};

  for(int ich=0; ich<5; ich++){ // 5 chambers per sector
    // Get real chamber number
    int chamber = (first_sector_chamber[sector-1] + ich);
    if(chamber > 18) chamber -= 18;

    int num_phi_addresses = std::pow(2,num_hs_bits + phi_addr_eta_bits);
    std::vector<int> chamber_phi_LUT; 
    chamber_phi_LUT.reserve(num_phi_addresses); // Need to fill whole addres space

    // Fill ME0 Phi LUT
    for(int addr=0; addr<num_phi_addresses; addr++){
      int phi_int = 0;

      // Get TP info from address
      int halfstrip = (addr >> phi_addr_eta_bits);
      int tp_partition_bitmask = (1 << phi_addr_eta_bits) - 1; // masks bits used as eta info (0x7)
      int tp_partition = (addr & tp_partition_bitmask) << (num_eta_bits - phi_addr_eta_bits); // Gets tp partition from address (lowest of possible ones)

      if(halfstrip < num_halfstrips){ // Rest of addresses are invalid
        // Get Chamber Partition
        int ch_partition = (tp_partition >> 1) + 1; // the actual partition in the chamber (1 indexed)
        ME0DetId detid(endcap, 0, chamber, ch_partition);
        const ME0Chamber* chamb = theME0Geometry_->chamber(detid);
        const ME0Layer* layer = chamb->layer(CSCConstants::KEY_ALCT_LAYER); // Eventually change this to ME0 key layer
        const ME0EtaPartition* roll = layer->etaPartition(ch_partition);

        // Get Global Point
        const float center_of_halfstrip = (0.5f * halfstrip) +  0.25; // Get center of halfstrip (0->1 is strip 1, so middle of hs 0 is 0.25)
        const LocalPoint& lp = roll->centreOfStrip(center_of_halfstrip);
        const GlobalPoint& gp = roll->surface().toGlobal(lp);

        float phi_degs = toolbox::rad_to_deg(gp.barePhi());
        phi_int = toolbox::calc_phi_int(phi_degs, sector);

        // Special ME0 case for chamber edge beyond sector 
        if(phi_int > 5040) // equivalent to ~84 degrees (this handle cases before and after sector edge bc of how colc_phi_int() works)
          phi_int = 0;
      }

      chamber_phi_LUT.push_back(phi_int);
    }


    int num_theta_addresses = std::pow(2,num_eta_bits);
    std::vector<int> chamber_theta_LUT; 
    chamber_theta_LUT.reserve(num_theta_addresses); // Need to fill whole addres space

    // Fill ME0 Theta LUT
    for(int addr=0; addr<num_theta_addresses; addr++){
      int tp_partition = addr;
      int ch_partition = (tp_partition >> 1) + 1; // the actual partition in the chamber (1 indexed)

      // Get Chamber Partition
      ME0DetId detid(endcap, 0, chamber, ch_partition);
      const ME0Chamber* chamb = theME0Geometry_->chamber(detid);
      const ME0Layer* layer = chamb->layer(CSCConstants::KEY_ALCT_LAYER); // Eventually change this to ME0 key layer
      const ME0EtaPartition* roll = layer->etaPartition(ch_partition);

      // Get Global Point
      const float strip = (0.25f * num_halfstrips);  // Gets halfstrip in ~middle of chamber (strip 192)
      const LocalPoint& lp = roll->centreOfStrip(strip);
      const GlobalPoint& gp = roll->surface().toGlobal(lp);

      // Get theta
      float theta_degs = toolbox::rad_to_deg(gp.theta());
      int theta_int = toolbox::calc_theta_int(theta_degs, endcap);
      chamber_theta_LUT.push_back(theta_int);
    }

    // Fill ME0 Phi LUT
    std::stringstream phi_filename; 
    phi_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/ME0_phi_ch" << ich << ".txt";
    std::ofstream phiFile(phi_filename.str());
    for (auto phi : chamber_phi_LUT) phiFile << phi << ", ";
    phiFile.close();

    // Fill ME0 Theta LUT
    std::stringstream theta_filename; 
    theta_filename << "PC_LUTs/endap" << endcap << "/sector" << sector << "/ME0_theta_ch" << ich << ".txt";
    std::ofstream thetaFile(theta_filename.str());
    for (auto theta : chamber_theta_LUT) thetaFile << theta << ", ";
    thetaFile.close();
  }

}



void MakePh2CoordLUT::generate_csc_LUTs(int endcap, int sector){

  int sector_idx = sector-1; // For start and end chs
  int start_chambers_10deg[6] = {2, 8, 14, 20, 26, 32};
  int start_chambers_20deg[6] = {1, 4, 7,  10, 13, 16};
  int end_chambers_10deg[6] = {9, 15, 21, 27, 33, 3}; // not inclusive
  int end_chambers_20deg[6] = {5, 8,  11, 14, 17, 2};

  // Using CSC Chamber Numbering
  for(int st=1; st<=4; st++){
    int max_ring = (st == 1)? 3 : 2; 
    for(int ring=1; ring<=max_ring; ring++){

      bool is_ten_deg = (st == 1 || ring == 2)? 1 : 0;
      int ch = (is_ten_deg)? start_chambers_10deg[sector_idx] : start_chambers_20deg[sector_idx];
      int end_ch = (is_ten_deg)? end_chambers_10deg[sector_idx] : end_chambers_20deg[sector_idx];
      int ich = 0; // chamber number within sector (0 is neighbor)

      while(ch != end_ch){
  
        if(st==1 && ring==1)
          generate_ME11_csc_LUT(endcap, sector, ch, ich);
        else
          generate_nonME11_csc_LUT(endcap, sector, st, ring, ch, ich);
       
        // Next chamber
        ich++;
        ch++;
        if((is_ten_deg && ch>36) || (!is_ten_deg && ch>18))
          ch = 1; // wrap to chamber 1
      } // ch
    } // ring
  } // st

}



void MakePh2CoordLUT::generate_ME11_csc_LUT(int endcap, int sector, int ch, int ich){
  
  int st = 1;

  // ME11b
  CSCDetId b_detid(endcap, st, 1, ch, CSCConstants::KEY_CLCT_LAYER); // Ring 1 - me11b
  const CSCChamber* b_chamb = theCSCGeometry_->chamber(b_detid);
  const CSCLayer* b_layer = b_chamb->layer(CSCConstants::KEY_ALCT_LAYER);
  const CSCLayerGeometry* b_layer_geom = b_layer->geometry();

  // ME11a
  CSCDetId a_detid(endcap, st, 4, ch, CSCConstants::KEY_CLCT_LAYER); // "Ring 4" - me11a
  const CSCChamber* a_chamb = theCSCGeometry_->chamber(a_detid);
  const CSCLayer* a_layer = a_chamb->layer(CSCConstants::KEY_ALCT_LAYER);
  const CSCLayerGeometry* a_layer_geom = a_layer->geometry();

  int me11a_strips = a_layer_geom->numberOfStrips();
  int me11b_strips = b_layer_geom->numberOfStrips();
  int num_wiregroups = b_layer_geom->numberOfWireGroups(); // Same for me11a and me11b

  int full_wg_bitwidth = 6; // 48 wiregroups
  int phi_addr_wg_bits = 3;
  int num_phi_addresses = std::pow(2, 10 + phi_addr_wg_bits); // 3 wg bits with 10 es bits
  int wg_interval = std::pow(2,(full_wg_bitwidth-phi_addr_wg_bits)); // ME11 - 8wiregroups per section

  std::vector<int> chamber_phi_LUT; 
  chamber_phi_LUT.reserve(num_phi_addresses); // Need to fill whole addres space

  // ----- Generate Phi LUT for Chamber -----------
  for(int addr=0; addr<num_phi_addresses; addr++){
    int phi_int = 0;

    int eighth_strip = (addr >> phi_addr_wg_bits); // remove wg bits to get eighth strip

    int wg_bit_mask = (1 << phi_addr_wg_bits) - 1; // masks bits used as wg (0x7)
    int fw_wg = (addr & wg_bit_mask) << (full_wg_bitwidth - phi_addr_wg_bits); // addr uses most significant wg bits, shift them into place (<<3) -> {0,8,16,24,32...}
    int wg = fw_wg + 1 + wg_interval/2; // Software wg 1 indexed, and add 4 to put wiregroup in middle of range
    

    // Note: in firmware, me11b strips come first (0-63), then me11a (64-...).
    // This is how the address space is organized
    // 
    // Wiregroup Sections for Phi calculation (sw indexing): 1-8, 9-16, 17-24, 25-32, 33-40, 41-48

    if(wg <= num_wiregroups){ 
      if( (eighth_strip/8) < me11b_strips){
        // In me11b address space
        GlobalPoint gp = get_csc_global_point(b_layer_geom, b_layer, wg, eighth_strip);
        float phi_degs = toolbox::rad_to_deg(gp.barePhi());
        if( fw_wg > 0) // wiregroups 1-8 don't overlap with me11b
          phi_int = toolbox::calc_phi_int(phi_degs, sector);
      }
      else if( (eighth_strip/8) < (me11b_strips + me11a_strips)){ // Check that 1/8th strip is in range
        // In me11a address space
        eighth_strip = eighth_strip - (me11b_strips * 8);
        GlobalPoint gp = get_csc_global_point(a_layer_geom, a_layer, wg, eighth_strip);
        float phi_degs = toolbox::rad_to_deg(gp.barePhi());
        if( fw_wg < 16) // wiregroups over 16 don't overlap with me11a
          phi_int = toolbox::calc_phi_int(phi_degs, sector);
      }
    }
      
    chamber_phi_LUT.push_back(phi_int);
  }


  // ----- Generate Theta LUT for Chamber -----------
  int th_addr_hs_bits = 6;
  int num_theta_addresses = std::pow(2, full_wg_bitwidth + th_addr_hs_bits); // 6 additional phi bits to deal with tilted wires
  std::vector<int> chamber_theta_LUT; 
  chamber_theta_LUT.reserve(num_theta_addresses); // Need to fill whole addres space

  for(int addr=0; addr<num_theta_addresses; addr++){ 
      int theta_int = 0; // Stays 0 for invalid values
      int wg = (addr >> th_addr_hs_bits) + 1; // remove phi bits to get wiregroup - 1 indexed so add 1

      int halfstrip_bitmask = (1 << th_addr_hs_bits) - 1; // masks bits used as hs (0x3f)
      int halfstrip = (addr & halfstrip_bitmask) << 2; // shift by 2 to get halfstrip bitwidth (8 bits)
      int eighth_strip = halfstrip << 2; // shift by 2 more to get 1/8th strip bitwidth of 10

      if(wg <= num_wiregroups){
        if( (eighth_strip/8) < me11b_strips){
          GlobalPoint gp = get_csc_global_point(b_layer_geom, b_layer, wg, eighth_strip);
          float theta_degs = toolbox::rad_to_deg(gp.theta());
          if(wg > 10)
            theta_int = toolbox::calc_theta_int(theta_degs, endcap);
        }
        else{
          GlobalPoint gp = get_csc_global_point(a_layer_geom, a_layer, wg, eighth_strip);
          float theta_degs = toolbox::rad_to_deg(gp.theta());
          if(wg < 17)
            theta_int = toolbox::calc_theta_int(theta_degs, endcap);
        }
      }
      chamber_theta_LUT.push_back(theta_int);
  }

  std::stringstream phi_filename; 
  phi_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/ME11_phi_ch" << ich << ".txt";
  std::ofstream phiFile(phi_filename.str());
  for (auto phi : chamber_phi_LUT) phiFile << phi << ", ";
  phiFile.close();

  std::stringstream theta_filename; 
  theta_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/ME11_theta_ch" << ich << ".txt";
  std::ofstream thetaFile(theta_filename.str());
  for (auto theta : chamber_theta_LUT) thetaFile << theta << ", ";
  thetaFile.close();
}



void MakePh2CoordLUT::generate_nonME11_csc_LUT(int endcap, int sector, int st, int ring, int ch, int ich){

  CSCDetId detid(endcap, st, ring, ch, CSCConstants::KEY_CLCT_LAYER);
  const CSCChamber* chamb = theCSCGeometry_->chamber(detid);
  const CSCLayer* layer = chamb->layer(CSCConstants::KEY_ALCT_LAYER);
  const CSCLayerGeometry* layer_geom = layer->geometry();

  int num_strips = layer_geom->numberOfStrips();
  int num_wiregroups = layer_geom->numberOfWireGroups();

  // Building Phi addresses with additional theta info
  int full_wg_bitwidth = std::ceil(std::log2(num_wiregroups)); // ring2=6bits, ring1=7bits, ME13=5bits
  int phi_addr_wg_bits = (ring == 1)? 2 : 1; // ring1=2wg bits, else=1wg bit
  int wg_interval = std::pow(2,(full_wg_bitwidth-phi_addr_wg_bits)); // ring1/2=32wgs, ME13 = 16wgs
  int max_wg_to_cover = std::pow(2,full_wg_bitwidth); // Need Addresses for all possible values in fw (even if wg is invalid)


  // ---------- Generate Phi LUT for chamber -----------
  int num_phi_addresses = std::pow(2,10 + phi_addr_wg_bits);
  std::vector<int> chamber_phi_LUT; 
  chamber_phi_LUT.reserve(num_phi_addresses); // Need to fill whole addres space

  for(int eighth_strip = 0; eighth_strip<std::pow(2,10); eighth_strip++){ // 0 is 1st 1/8th strip of first strip (strip 1 in csc terms)
    for(int wg=(wg_interval/2); wg<=max_wg_to_cover; wg+=wg_interval){ // use wg in middle of interval range
      int phi_int = 0; // invalid value

      if(wg <= num_wiregroups && eighth_strip/8 < num_strips){
        GlobalPoint gp = get_csc_global_point(layer_geom, layer, wg, eighth_strip);
        float phi_degs = toolbox::rad_to_deg(gp.barePhi());
        phi_int   = toolbox::calc_phi_int(phi_degs, sector);
      }

      chamber_phi_LUT.push_back(phi_int);
    }
  }


  // ----- Generate Theta LUT for Chamber -----------
  int num_theta_addresses = std::pow(2, full_wg_bitwidth); 
  std::vector<int> chamber_theta_LUT; 
  chamber_theta_LUT.reserve(num_theta_addresses); // Need to fill whole addres space

  for(int addr=0; addr<num_theta_addresses; addr++){ 
      int theta_int = 0;
      int wg = addr + 1; // wiregroups start at 1 for csc geometry functions
      int eighth_strip = (num_strips * 8) / 2; // use phi midpoint in chamber

      if(wg <= num_wiregroups){
        GlobalPoint gp = get_csc_global_point(layer_geom, layer, wg, eighth_strip);
        float theta_degs = toolbox::rad_to_deg(gp.theta());
        theta_int = toolbox::calc_theta_int(theta_degs, endcap);
      }
      
      chamber_theta_LUT.push_back(theta_int);
  }

  std::stringstream phi_filename; 
  phi_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/ME" << st << ring << "_phi_ch" << ich << ".txt";
  std::ofstream phiFile(phi_filename.str());
  for (auto phi : chamber_phi_LUT) phiFile << phi << ", ";
  phiFile.close();

  std::stringstream theta_filename; 
  theta_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/ME" << st << ring << "_theta_ch" << ich << ".txt";
  std::ofstream thetaFile(theta_filename.str());
  for (auto theta : chamber_theta_LUT) thetaFile << theta << ", ";
  thetaFile.close();
}



GlobalPoint MakePh2CoordLUT::get_csc_global_point(const CSCLayerGeometry* layer_geom, const CSCLayer* layer, const int wiregroup, const int eighth_strip){
  // Convert to float strip and integer wire number
  float middle_of_eighth_strip = (eighth_strip + 0.5) / 8;
  int   middle_wire_of_group = int(layer_geom->middleWireOfGroup(wiregroup)); // For even number of wires in group (causing floating point wire number), this number is truncated

  const LocalPoint& lp =  layer_geom->intersectionOfStripAndWire( middle_of_eighth_strip, middle_wire_of_group); // this function uses float strip and integer wire
  const GlobalPoint& gp = layer->surface().toGlobal(lp);

  return gp;
}



void MakePh2CoordLUT::quick_LUT_verification(int endcap, int sector){
  int ph_hard_init_10deg[7] = {38, 75, 113, 150, 188, 225, 263};
  int ph_hard_init_20deg[4] = {0,  75, 150, 225};
  int ph_hard_init_me0[5] =   {0,  56, 131, 206, 281};

  int cover_10deg = 53;
  int cover_20deg = 90;
  int cover_me0 = 90;

  // Check ME11
  for(int ich=0; ich<7; ich++){
    std::stringstream phi_filename;
    phi_filename << "PC_LUTs/endcap" << endcap << "/sector" << sector << "/ME11_phi_ch" << ich << ".txt";
    std::ifstream phiFile(phi_filename.str());  
    
    int val, col_val;
    char comma;
    int addr = 0;
    while(phiFile >> val >> comma){
      col_val = (val >> 4);
      if( val != 0 && (col_val < ph_hard_init_10deg[ich] || (col_val > (ph_hard_init_10deg[ich] + cover_10deg))))
        std::cout << "ERROR - VALUE IS OUT OF RANGE - ADDR: " << addr << "\tValue: " << val << " - " << col_val << std::endl;
      addr++;
    }
  
    std::cout << "Final Address: " << addr << std::endl;
  }

}






// DEFINE THIS AS A PLUG-IN
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakePh2CoordLUT);
