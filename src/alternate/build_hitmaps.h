#ifndef __PROMPT_TRKBUILD_ZONING_H__
#define __PROMPT_TRKBUILD_ZONING_H__

//#include "../emtf_hlslib/layer_helpers.h"
#include <iostream>
#include "common.h"
#include "types.h"

namespace prompt_trkbuild {

void build_hitmaps( emtf_seg_t emtf_segs[num_sites][max_ch_site][max_seg_ch], hitmap_t hitmaps[num_zones][num_hitmap_rows] ){

    raw_hitmap_t raw_hitmap[num_zones][num_hitmap_rows];

    // Set all bits to 0
    for(int zone=0; zone<num_zones; zone++)
        for(int row=0; row<num_hitmap_rows; row++)
            raw_hitmap[zone][row] = 0;

    // go through each zone and row and create hitmap for that row
    for(int zone=0; zone<num_zones; zone++)
	    for(int row=0; row<num_hitmap_rows; row++)

            // Build hitmap by looping through all segments in a row and setting column corresponding to segment phi = 1
            for(int site=0; site<max_sites_in_row; site++){
                for(int ch=0; ch<max_ch_site; ch++)
                    for(int seg=0; seg<max_seg_ch; seg++){

                        int curr_site = sites_by_row[zone][row][site];
                        if(curr_site != no_site and emtf_segs[curr_site][ch][seg].seg_valid and emtf_segs[curr_site][ch][seg].seg_zones[2-zone]){
                            unsigned int col = (emtf_segs[curr_site][ch][seg].emtf_phi >> 4);
                            raw_hitmap[zone][row][col] = (raw_hitmap[zone][row][col] | 1);
                        }
                    }
            }


    // Set output to last 288 bits of raw hitmap
    for(int zone=0; zone<num_zones; zone++)
        for(int row=0; row<num_hitmap_rows; row++)
            hitmaps[zone][row] = raw_hitmap[zone][row].range(314,27);

}


} // namespace


#endif
