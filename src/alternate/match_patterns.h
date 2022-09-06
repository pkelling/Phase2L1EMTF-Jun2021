#ifndef __PROMPT_TRKBUILD_PATT_MATCH_H__
#define __PROMPT_TRKBUILD_PATT_MATCH_H__

//#include "../emtf_hlslib/layer_helpers.h"
#include <iostream>
#include "common.h"
#include "types.h"
#include "patterns.h"

namespace prompt_trkbuild {


void match_patterns(hitmap_t hitmaps[num_zones][num_hitmap_rows], patt_match_t patt_matches[num_zones][num_hitmap_cols] ){

    ap_uint<num_hitmap_rows> activations[num_zones][num_patterns][num_hitmap_cols];	// one for each column

    for(int zone=0; zone < num_zones; zone++){
        for(int row=0; row < num_hitmap_rows; row++){
            for(int patt=0; patt < num_patterns; patt++){

                // Pad the hitmap with zeros to cover cases where pattern range is out of range (288 bits + max padding on each side)
                ap_uint<num_hitmap_cols+max_padding*2> padded_hitmap_row = (ap_uint<max_padding>(0), hitmaps[zone][row], ap_uint<max_padding>(0));

                // Slide the patterns across the hitmap and check for 'activations' (matches)
                for(int col=0; col<num_hitmap_cols; col++){
                    int col_start = col + col_start_table[zone][patt][row];
                    int col_stop = col + col_stop_table[zone][patt][row];

                    // Bitwise or for hitmap cols in pattern range
                    activations[zone][patt][col][row] = padded_hitmap_row.range(col_stop+1,col_start).or_reduce();
                }
		    }
	    }
    }

    // Get activation qualities
    ap_uint<6> activation_quality[num_zones][num_patterns];
    ap_uint<6> max_quality = 0;
    ap_uint<3> max_patt = 0;

    for(int zone=0; zone < num_zones; zone++){
        for(int col=0; col<288; col++){
            for(int patt=0; patt < 7; patt++){
                ap_uint<num_hitmap_rows> col_activation = activations[zone][patt][col];
                activation_quality[zone][patt] = pattern_activations[zone][col_activation];
            }

            // Sort Activations
            max_quality = 0;
            max_patt = 0;

            for(int patt=0; patt<7; patt++){
                // find highest quality for this column (priority goes to lower pattern number)
                if(activation_quality[zone][patt] > max_quality){
                    max_quality = activation_quality[zone][patt];
                    max_patt = patt;
                }
            }

            patt_matches[zone][col].patt = max_patt;
            patt_matches[zone][col].qual = max_quality;
        }
    }

}


} // namespace


#endif
