#ifndef __PROMPT_TRKBUILD_INPUT_CONV_H__
#define __PROMPT_TRKBUILD_INPUT_CONV_H__

//#include "../emtf_hlslib/layer_helpers.h"
#include <iostream>
#include "common.h"
#include "types.h"

namespace prompt_trkbuild {


void inline convert_input(  emtf_phi_t phi[230], 
                            emtf_theta_t th1[230], 
                            emtf_theta_t th2[230], 
                            emtf_bend_t bend[230], 
                            emtf_qual_t qual[230], 
                            emtf_time_t seg_times[230],
                            seg_zones_t zones[230],
                            ap_uint<3> 	tzones[230],
                            seg_bx_t 	bx[230],
                            seg_valid_t valid[230],
                            emtf_seg_t emtf_segs[num_sites][max_ch_site][max_seg_ch]
                        )
{
    
    
    int ch_ids[num_sites][max_ch_site] = {
        
        { 108, 109, 110, 111, 112, 113, 114},   // ME0      // This will change (only 5 chambers)
        { 54,  55,  56,  63,  64,  65,  99},    // GE11
        { 0,   1,   2,   9,  10,  11,  45},     // ME11
        { 3,   4,   5,  12,  13,  14,  46},     // ME12
        { 57,  58,  59,  66,  67,  68, 100},    // RE12
        { 72,  73,  74, 102, -1, -1, -1},       // GE21
        { 75,  76,  77,  78,  79,  80, 103},    // RE22
        { 18,  19,  20,  48, -1, -1, -1},       // ME21
        { 21,  22,  23,  24,  25,  26,  49},    // ME22
        { 27,  28,  29,  50, -1, -1, -1},       // ME31
        { 30,  31,  32,  33,  34,  35,  51},    // ME32
        { 81,  82,  83, 104, -1, -1, -1},       // RE31
        { 84,  85,  86,  87,  88,  89, 105},    // RE32
        { 36,  37,  38,  52, -1, -1, -1},       // ME41
        { 39,  40,  41,  42,  43,  44,  53},    // ME42
        { 90,  91,  92, 106, -1, -1, -1},       // RE41
        { 93,  94,  95,  96,  97,  98, 107}     // RE42
    };

    for(int site=0; site<num_sites; site++)
        for(int ch=0; ch<max_ch_site; ch++)
            for(int seg=0; seg<max_seg_ch; seg++){
                int ch_id = ch_ids[site][ch];

                if(ch_id != -1 && seg<2){
                    emtf_segs[site][ch][seg].emtf_phi    = phi[2*ch_id + seg];
                    emtf_segs[site][ch][seg].emtf_theta1 = th1[2*ch_id + seg];
                    emtf_segs[site][ch][seg].emtf_theta2 = th2[2*ch_id + seg];
                    emtf_segs[site][ch][seg].emtf_bend   = bend[2*ch_id + seg];
                    emtf_segs[site][ch][seg].emtf_qual   = qual[2*ch_id + seg];
                    emtf_segs[site][ch][seg].emtf_time   = seg_times[2*ch_id + seg];
                    emtf_segs[site][ch][seg].seg_zones  = zones[2*ch_id + seg];
                    emtf_segs[site][ch][seg].seg_bx     = bx[2*ch_id + seg];
                    emtf_segs[site][ch][seg].seg_valid  = valid[2*ch_id + seg];

                    if( tzones[2*ch_id + seg][2] != 1)
                        emtf_segs[site][ch][seg].seg_valid  = 0;
                }
                else{
                    emtf_segs[site][ch][seg].emtf_phi    = 0;
                    emtf_segs[site][ch][seg].emtf_theta1 = 0;
                    emtf_segs[site][ch][seg].emtf_theta2 = 0;
                    emtf_segs[site][ch][seg].emtf_bend   = 0;
                    emtf_segs[site][ch][seg].emtf_qual   = 0;
                    emtf_segs[site][ch][seg].emtf_time   = 0;
                    emtf_segs[site][ch][seg].seg_zones  = 0;
                    emtf_segs[site][ch][seg].seg_bx     = 0;
                    emtf_segs[site][ch][seg].seg_valid  = 0;
                }
            }
}


} // namespace


#endif
