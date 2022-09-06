#ifndef __PROMPT_TRKBUILD_FIND_BEST_TRACK_H__
#define __PROMPT_TRKBUILD_FIND_BEST_TRACK_H__

//#include "../emtf_hlslib/layer_helpers.h"
#include <iostream>
#include "common.h"
#include "types.h"

namespace prompt_trkbuild {

template <typename T>
void sort_four_op(const T in[4], T out[4]) {
    // Sorts four unsorted inputs (not stable) - input struct must have .qual

    // Stage 1: compare-swap if (wire[i] < wire[j]) swap(wire[j], wire[i])
    const T tmp_1_0 = (in[0].qual < in[1].qual) ? in[1] : in[0];
    const T tmp_1_1 = (in[0].qual < in[1].qual) ? in[0] : in[1];
    const T tmp_1_2 = (in[2].qual < in[3].qual) ? in[3] : in[2];
    const T tmp_1_3 = (in[2].qual < in[3].qual) ? in[2] : in[3];

    // Stage 2
    const T tmp_2_0 = (tmp_1_0.qual < tmp_1_2.qual) ? tmp_1_2 : tmp_1_0;
    const T tmp_2_1 = (tmp_1_1.qual < tmp_1_3.qual) ? tmp_1_3 : tmp_1_1;
    const T tmp_2_2 = (tmp_1_0.qual < tmp_1_2.qual) ? tmp_1_0 : tmp_1_2;
    const T tmp_2_3 = (tmp_1_1.qual < tmp_1_3.qual) ? tmp_1_1 : tmp_1_3;

    // Stage 3
    const T tmp_3_1 = (tmp_2_1.qual < tmp_2_2.qual) ? tmp_2_2 : tmp_2_1;
    const T tmp_3_2 = (tmp_2_1.qual < tmp_2_2.qual) ? tmp_2_1 : tmp_2_2;

    // Output
    out[0] = tmp_2_0;
    out[1] = tmp_3_1;
    out[2] = tmp_3_2;
    out[3] = tmp_2_3;
}


template <typename T>
void merge_eight_op(const T in[8], T out[4]){
    // Bitonic Parrallel sort, first four and last four must be pre-sorted - input struct must have .qual
    // This is not "stable" (will not keep order of same values)

    // Stage 1: compare-swap if (wire[i] < wire[j]) swap(wire[j], wire[i])
    const T tmp_1_0 = (in[0].qual < in[4].qual) ? in[4] : in[0];
    const T tmp_1_1 = (in[1].qual < in[5].qual) ? in[5] : in[1];
    const T tmp_1_2 = (in[2].qual < in[6].qual) ? in[6] : in[2];
    const T tmp_1_3 = (in[3].qual < in[7].qual) ? in[7] : in[3];
    const T tmp_1_4 = (in[0].qual < in[4].qual) ? in[0] : in[4];
    const T tmp_1_5 = (in[1].qual < in[5].qual) ? in[1] : in[5];

    // Stage 2
    const T tmp_2_2 = (tmp_1_2.qual < tmp_1_4.qual) ? tmp_1_4 : tmp_1_2;
    const T tmp_2_3 = (tmp_1_3.qual < tmp_1_5.qual) ? tmp_1_5 : tmp_1_3;
    const T tmp_2_4 = (tmp_1_2.qual < tmp_1_4.qual) ? tmp_1_2 : tmp_1_4;

    // Stage 3
    const T tmp_3_1 = (tmp_1_1.qual < tmp_2_2.qual) ? tmp_2_2 : tmp_1_1;
    const T tmp_3_2 = (tmp_1_1.qual < tmp_2_2.qual) ? tmp_1_1 : tmp_2_2;
    const T tmp_3_3 = (tmp_2_3.qual < tmp_2_4.qual) ? tmp_2_4 : tmp_2_3;

    // Output
    out[0] = tmp_1_0;
    out[1] = tmp_3_1;
    out[2] = tmp_3_2;
    out[3] = tmp_3_3;
    
}


// Input: (3bit patt, 6bit qual)
// Output: (2bit zone, 8bit col, 3bit patt, 6bit qual)

void find_best_tracks(patt_match_t patt_matches[num_zones][num_hitmap_cols], best_trk_t best_trks[num_tracks]){


    // Intermediate Output  - Order: [ zone 0 trks, zone 1 trks, zone 2 trks]
    best_trk_t best_zone_trks[num_zones * num_tracks];


    // First, find the best tracks from each zone individually
    for(int zone=0; zone<num_zones; zone++){

        // Step 1: suppress if not local maximum    if(qcol <= qbefore || qcol < qafter)   
        best_trk_t suppressed_qualities[num_hitmap_cols];
        for(int col=0; col<num_hitmap_cols; col++){

            if (    (col != (num_hitmap_cols-1) && patt_matches[zone][col].qual < patt_matches[zone][col+1].qual) or
                    (col != 0                   && patt_matches[zone][col].qual <= patt_matches[zone][col-1].qual)){
                suppressed_qualities[col].qual = 0;
            }
            else{
                suppressed_qualities[col].col = col;
                suppressed_qualities[col].qual = patt_matches[zone][col].qual;
                suppressed_qualities[col].patt = patt_matches[zone][col].patt;
            }
        }
        

        // Pick best of every neighbor pair
        best_trk_t muxed_qualities[num_hitmap_cols/2]; 
        for(int i=0; i<(num_hitmap_cols/2); i++)
            muxed_qualities[i] = (suppressed_qualities[2*i].qual > 0) ? suppressed_qualities[2*i] : suppressed_qualities[2*i+1];


        // Setup for sorting - sort every four values of first stage (zero all other stage qualities)
        best_trk_t data_stages[7][num_hitmap_cols];
        for(int stage=0; stage<7; stage++){
            for(int mcol=0; mcol<(num_hitmap_cols/2); mcol+=4){ // muxed column
                if(stage == 0){
                    // Sorts every four values for first stage
                    sort_four_op( &muxed_qualities[mcol], &data_stages[0][mcol]);
                }
                else{
                    data_stages[stage][mcol].qual = 0;
                    data_stages[stage][mcol+1].qual = 0;
                    data_stages[stage][mcol+2].qual = 0;
                    data_stages[stage][mcol+3].qual = 0;
                }
		    }
        }


        // Skipping nodes is done to make it equivalent to Jia Fu's 
        int nodes[6] = {18, 9, 5, 3, 2, 1};	// including possible extra stage with no pair
        int node_to_skip[6] = {1000, 1000, 1, 1, 1, 1000};

        // Main Sorting for each zone
        for(int stage=0; stage<6; stage++){
            for(int node=0; node<nodes[stage]; node+=1){
                // Bitonic Sort 8 into best 4 (each set of 4 is already sorted)

                int mcol = 8 * node;

                if(node == node_to_skip[stage]){
                    // Output
                    data_stages[stage+1][mcol/2] =  	data_stages[stage][mcol];
                    data_stages[stage+1][mcol/2+1] =  data_stages[stage][mcol+1];
                    data_stages[stage+1][mcol/2+2] =  data_stages[stage][mcol+2];
                    data_stages[stage+1][mcol/2+3] =  data_stages[stage][mcol+3];
                }
                else{
                    if(node > node_to_skip[stage])
                        mcol-=4;

                    // Select top 4 from every 8
                    merge_eight_op( &data_stages[stage][mcol], &data_stages[stage+1][4*node]);
                }
            }
        }

        // Collect intermediate outputs
        for(int trk=0; trk<num_tracks; trk++){
            best_zone_trks[zone * num_tracks + trk].qual = data_stages[6][trk].qual;
            best_zone_trks[zone * num_tracks + trk].patt = data_stages[6][trk].patt;
            best_zone_trks[zone * num_tracks + trk].col = data_stages[6][trk].col;
            best_zone_trks[zone * num_tracks + trk].zone = zone;
        }

    } // zone

    best_trk_t data_stages2[3][16];
    for(int stage=0; stage<2; stage++){
	    for(int trk=0; trk<16; trk++){
		    if(stage == 0 and trk<12)
                data_stages2[stage][trk] = best_zone_trks[trk]; // z0, z1, then z2
		    else
			    data_stages2[stage][trk].qual = 0;
	    }
    }
    
    /************* Zone Merging (sort between top tracks from each zone) ******************/
    int nodes[2] = {2, 1};	// including possible extra stage with no pair
    for(int stage=0; stage<2; stage++){
        for(int node=0; node<nodes[stage]; node+=1){
            int mcol = 8 * node;
            merge_eight_op( &data_stages2[stage][mcol], &data_stages2[stage+1][4*node]);
        }
    }

    // Collect outputs
    for(int trk=0; trk<num_tracks; trk++)
        best_trks[trk] = data_stages2[2][trk];

}


} // namespace


#endif
