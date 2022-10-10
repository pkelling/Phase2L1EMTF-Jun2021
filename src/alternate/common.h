#ifndef __PROMPT_TRKBUILD_COMMON_H__
#define __PROMPT_TRKBUILD_COMMON_H__

#include <string>

namespace prompt_trkbuild {

// EMTF specific
constexpr int num_zones = 3;      // per sector

constexpr int num_sites = 17;       // per track
constexpr int max_ch_site = 7;

enum Sites {ME0=0, GE11, ME11, ME12, RE12, GE21, RE22, ME21, ME22, ME31, ME32, RE31, RE32, ME41, ME42, RE41, RE42=16 };

// Only for printing
const std::string site_strings[num_sites] = {
            "ME0", "GE11", "ME11", "ME12", "RE12", "GE21", "RE22", "ME21", "ME22", "ME31", "ME32", "RE31", "RE32", "ME41", "ME42", "RE41", "RE42" 
};


// Segments Per Chamber
constexpr int max_seg_ch = 8;
constexpr int seg_ch_csc = 2;   
constexpr int seg_ch_rpc12 = 2; // rpc stations 1 and 2
constexpr int seg_ch_rpc34 = 4; // rpc stations 3 and 4 (combo chambers)
constexpr int seg_ch_ge11 = 4;   
constexpr int seg_ch_ge21 = 4;   
constexpr int seg_ch_me0 = 6;


constexpr int num_tracks = 4;       // per sector
constexpr int num_patterns = 7;     // per zone
constexpr int num_emtf_features = 40;    // per track

constexpr int num_hitmap_rows = 8;
constexpr int num_raw_hitmap_cols = 315;
constexpr int num_hitmap_cols = 288;


const int max_sites_in_row = 4;
constexpr int max_num_layers = 12;
const int max_chambers_in_layer = 11; // Layer for ring 1 and ring 2


// Chamber Types
const int num_site_types = 3;
enum Angle_Types { ten_deg=0, twenty_deg=1, twenty_me0=2 };

const int unsigned type_by_site[num_sites] = { 
                                            ten_deg, ten_deg, ten_deg, // twenty_me0, ten_deg, ten_deg,               // Station 0 
                                            ten_deg, ten_deg,                           // Station 1
                                            twenty_deg, ten_deg, twenty_deg, ten_deg,   // Station 2
                                            twenty_deg, ten_deg, twenty_deg, ten_deg,   // Station 3
                                            twenty_deg, ten_deg, twenty_deg, ten_deg    // Staiton 4
};

// Hitmap Information

const int no_site = -1;
const int sites_by_row[num_zones][num_hitmap_rows][max_sites_in_row] = 
{ 
    // Zone 0
    {
        {ME0,  no_site, no_site, no_site}, 
        {GE11, no_site, no_site, no_site}, 
        {ME11, no_site, no_site, no_site}, 
        {GE21, no_site, no_site, no_site}, 
        {ME21, no_site, no_site, no_site},  
        {ME31, no_site, no_site, no_site},  
        {RE31, no_site, no_site, no_site}, 
        {ME41, RE41,    no_site, no_site}
    },
    
    // Zone 1
    { 
        {GE11, no_site, no_site, no_site}, 
        {ME11, no_site, no_site, no_site},  
        {ME12, RE12,    no_site, no_site}, 
        {GE21, no_site, no_site, no_site},  
        {ME21, no_site, no_site, no_site},  
        {ME31, ME32,    no_site, no_site}, 
        {RE31, RE32,    no_site, no_site}, 
        {ME41, ME42,    RE41,    RE42}
    },
    
    // Zone 2
    { 
        {ME12, no_site, no_site, no_site}, 
        {RE12, no_site, no_site, no_site}, 
        {RE22, no_site, no_site, no_site}, 
        {ME22, no_site, no_site, no_site},  
        {ME32, no_site, no_site, no_site}, 
        {RE32, no_site, no_site, no_site}, 
        {ME42, no_site, no_site, no_site}, 
        {RE42, no_site, no_site, no_site}
    }                                                                                                                                               
}; 


// Some sites are combined into 1 layer for the NN (ordered based on Jia Fus NN)
const int map_layers_to_sites[max_num_layers][2] = { 
                                    {ME11, no_site}, 
                                    {ME12, no_site},
                                    {ME21, ME22},
                                    {ME31, ME32},
                                    {ME41, ME42},
                                    {RE12, no_site},
                                    {RE22, no_site},
                                    {RE31, RE32},
                                    {RE41, RE42},
                                    {GE11, no_site},
                                    {GE21, no_site},
                                    {ME0, no_site}
};


const int map_layers_to_zone_rows[num_zones][max_num_layers] = {
                { 2, -1, 4, 5, 7, -1, -1, 6, 7,  1,  3,  0},
                { 1,  2, 4, 5, 7,  2, -1, 6, 7,  0,  3, -1},
                { -1, 0, 3, 4, 6,  1,  2, 5, 7, -1, -1, -1}
};


}  // namespace 




#endif
