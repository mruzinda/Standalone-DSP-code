/* flag_xrfi_thread.c
 * 
 * Routine to form beams and mitigate RFI from received packets
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "xrfi.h"
#include "hashpipe.h"
#include "flag_databuf.h"

// Create thread status buffer
static hashpipe_status_t * st_p;

typedef enum {
    ACQUIRE,
    CLEANUP
} state;


// Run method for the thread
// It is meant to do the following:
//     (1) Initialize status buffer
//     (2) Start main loop
//         (2a) Wait for input buffer block to be filled
//         (2b) Print out some data in the block
static void * run(hashpipe_thread_args_t * args) {
    // Local aliases to shorten access to args fields
    flag_input_databuf_t * db_in = (flag_input_databuf_t *)args->ibuf;
    flag_gpu_beamformer_output_databuf_t * db_out = (flag_gpu_beamformer_output_databuf_t *)args->obuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    st_p = &st; // allow global (this source file) access to the status buffer

    int rv;
    int curblock_in = 0;
    int curblock_out = 0;
    uint64_t start_mcnt = 0;
    uint64_t last_mcnt = Nm - 1;

#if XRFI_TIME
    float time_taken = 0;
    float time_taken1 = 0;
    float time_taken2 = 0;
    float corr_time = 0;
    float proj_time = 0;
    float beam_time = 0;
#endif

    // Initialize XRFI filter
    printf("Initalizing XRFI filter...\n\r");
    init_xrfi();

    // Initialize beamformer
    printf("Initializing beamformer...\n\r");
    init_beamformer_xrfi();

    // Update weights
    // TODO: allow update of weights during runtime
    printf("RTB: Initializing beamformer weights...\n");
    char weightdir[65];
    char weight_file[17];
    hashpipe_status_lock_safe(&st);
    hgets(st.buf,"WEIGHTD", 65, weightdir);
    hashpipe_status_unlock_safe(&st);

    hashpipe_status_lock_safe(&st);
    hgets(st.buf,"BWEIFILE",16,weight_file);
    hashpipe_status_unlock_safe(&st);

    char weight_flag[9];
    hashpipe_status_lock_safe(&st);
    hgets(st.buf, "WFLAG", 8, weight_flag);
    hashpipe_status_unlock_safe(&st);
  
    char w_dir[70];
    //sprintf(w_dir, "%s/dummy_w_A.bin", weightdir);
    printf("RTBF: Weight flag: %s.\n", weight_flag);
    if (strncmp(weight_flag, "1",1) == 0){
        printf("RTBF: Updated weights.\n");
        sprintf(w_dir, "%s/%s", weightdir, weight_file);
    }
    else {
        printf("RTBF: Unity weights used.\n");
        sprintf(w_dir, "%s/dummy_w_A.bin", weightdir);
    }
    printf("BF: Weight file name: %s\n", w_dir);

    // update_weights("./weights.in");
    update_weights(w_dir);
    // Put metadata into status shared memory
    float offsets[BX_BEAM];
    char cal_filename[BX_META];
    char algorithm[BX_META];
    char weight_filename[BX_META];
    //long long unsigned int bf_xid = 0;
    int act_xid;

    bf_get_offsets(offsets);
    bf_get_cal_filename(cal_filename);
    bf_get_algorithm(algorithm);
    bf_get_weight_filename(weight_filename);
    //long long unsigned int bf_xid = bf_get_xid();

    int i;
    hashpipe_status_lock_safe(&st);
    for (i = 0; i < BX_BEAM/2; i++) {
        char keyword1[9];
        snprintf(keyword1,8,"ELOFF%d",i);
        hputr4(st.buf, keyword1, offsets[2*i]);
        char keyword2[9];
        snprintf(keyword2,8,"AZOFF%d",i);
        hputr4(st.buf, keyword2, offsets[2*i+1]);
    }
    hputs(st.buf, "BCALFILE", cal_filename);
    hputs(st.buf, "BALGORIT", algorithm);
    hputs(st.buf, "BWFILE", weight_filename);
    hgeti4(st.buf, "XID", &act_xid);
    hashpipe_status_unlock_safe(&st);
    printf("RTBF: Finished weight initialization at beginning of thread...\n");

    state cur_state = ACQUIRE;
    state next_state = ACQUIRE;

    char netstat[17];

    // Indicate in shared memory buffer that this thread is ready to start
    hashpipe_status_lock_safe(&st);
    hputi4(st.buf, "RBFREADY", 1);
    hashpipe_status_unlock_safe(&st);
    
    printf("Using host arrays allocated in pinned memory\n\r");

    for (i = 0; i < N_INPUT_BLOCKS; i++){
        ///////////////////////////////////////////////////////////////////////////////
        //>>>>   Register host array in pinned memory <<<<
        ///////////////////////////////////////////////////////////////////////////////
        input_data((signed char *)&db_in->block[i].data);
    }
    for (i = 0; i < N_GPU_OUT_BLOCKS; i++){
        ///////////////////////////////////////////////////////////////////////////////
        //>>>>   Register host array in pinned memory <<<<
        ///////////////////////////////////////////////////////////////////////////////
        output_data((float *)&db_out->block[i].data);
    }

    ////////////////////////////////////////////////////////////////
    //>>>>    Allocating output host array in pinned memory   <<<<
    ////////////////////////////////////////////////////////////////
    float * c_data_out;
    c_data_out = (float *)calloc(2*BX_ELE*BX_ELE*BX_BIN, sizeof(float));
    corr_output_data(c_data_out);

    // Main loop for thread
    while (run_threads()) {
        if(cur_state == ACQUIRE){
            next_state = ACQUIRE;
	        // Wait for input buffer block to be filled
            while ((rv=flag_input_databuf_wait_filled(db_in, curblock_in)) != HASHPIPE_OK && run_threads()) {
                if (rv==HASHPIPE_TIMEOUT) {
                    int cleanA;
                    hashpipe_status_lock_safe(&st);
                    hgetl(st.buf, "CLEANA", &cleanA);
                    hgets(st.buf, "NETSTAT", 16, netstat);
                    hashpipe_status_unlock_safe(&st);
                    if (cleanA == 0 && strcmp(netstat, "CLEANUP") == 0) {
                        next_state = CLEANUP;
                        break;
                    }
                }
                else {
                    hashpipe_error(__FUNCTION__, "error waiting for filled databuf block");
                    pthread_exit(NULL);
                    break;
                }
            }
            if (!run_threads()) break;

            // If CLEANUP, don't continue processing
            if (next_state != CLEANUP) {

                if (DEBUG) {
                    // Print out the header information for this block 
                    flag_input_header_t tmp_header;
                    memcpy(&tmp_header, &db_in->block[curblock_in].header, sizeof(flag_input_header_t));
                    hashpipe_status_lock_safe(&st);
                    hputi4(st.buf, "BEAMMCNT", tmp_header.mcnt_start);
                    hashpipe_status_unlock_safe(&st);
                }

                // Wait for output block to become free
                while ((rv=flag_gpu_beamformer_output_databuf_wait_free(db_out, curblock_out)) != HASHPIPE_OK) {
                    if (rv==HASHPIPE_TIMEOUT) {
                        continue;
                    } else {
                        hashpipe_error(__FUNCTION__, "error waiting for free databuf");
                        fprintf(stderr, "rv = %d\n", rv);
                        pthread_exit(NULL);
                        break;
                    }
                }

                ///////////////////////////////////////////////////////////
                //>>>>    Run correlator   <<<<
                ///////////////////////////////////////////////////////////
#if VERBOSE
                printf("Running correlator \n\r");
#endif

#if XRFI_TIME
		//struct timeval tval_before, tval_after, tval_result;
		//gettimeofday(&tval_before, NULL);
		struct timespec tval_before, tval_after;
		clock_gettime(CLOCK_MONOTONIC, &tval_before);
#endif
		/////////////////////////////////////////////////////////////
		run_correlator_xrfi((signed char *)&db_in->block[curblock_in].data, c_data_out);
		//run_correlator_xrfi(c_data_out);
		/////////////////////////////////////////////////////////////

#if XRFI_TIME
		//gettimeofday(&tval_after, NULL);
		//timersub(&tval_after, &tval_before, &tval_result);
		//corr_time = (float)tval_result.tv_usec/1000;
		clock_gettime(CLOCK_MONOTONIC, &tval_after);
		time_taken = (float)(tval_after.tv_sec - tval_before.tv_sec)*1e6; // Time in seconds since epoch
		time_taken = time_taken + (float)(tval_after.tv_nsec - tval_before.tv_nsec)*1e-6; // Time in nanoseconds since 'tv_sec - start and end'
		corr_time = time_taken;
                //free(c_data_out);
#endif
                ////////////////////////////////////////////////////////////

                ///////////////////////////////////////////////////////////
                //>>>>    Update Weights for RFI mitigation   <<<<
                ///////////////////////////////////////////////////////////
#if VERBOSE
		printf("Computing Subspace\n\r");
#endif

#if XRFI_TIME
		//struct timeval tval_before1, tval_after1, tval_result1;
		//gettimeofday(&tval_before1, NULL);
		struct timespec tval_before1, tval_after1;
		clock_gettime(CLOCK_MONOTONIC, &tval_before1);
#endif

		xrfiComputeBatchedProjectionMatrix();

#if VERBOSE
		printf("Subspace computed\n\r");
#endif

#if XRFI_TIME
		//gettimeofday(&tval_after1, NULL);
		//timersub(&tval_after1, &tval_before1, &tval_result1);
		//proj_time = (float)tval_result1.tv_usec/1000;
		clock_gettime(CLOCK_MONOTONIC, &tval_after1);
		time_taken1 = (float)(tval_after1.tv_sec - tval_before1.tv_sec)*1e6; // Time in seconds since epoch
		time_taken1 = time_taken1 + (float)(tval_after1.tv_nsec - tval_before1.tv_nsec)*1e-6; // Time in nanoseconds since 'tv_sec - start and end'
		proj_time = time_taken1;
#endif
		////////////////////////////////////////////////////////////

                ///////////////////////////////////////////////////////////
                //>>>>    Run Beamformer   <<<<
                ///////////////////////////////////////////////////////////
#if XRFI_TIME
		//struct timeval tval_before2, tval_after2, tval_result2;
		//gettimeofday(&tval_before2, NULL);
		struct timespec tval_before2, tval_after2;
		clock_gettime(CLOCK_MONOTONIC, &tval_before2);
#endif

                run_beamformer_xrfi((signed char *)&db_in->block[curblock_in].data, (float *)&db_out->block[curblock_out].data);
                //run_beamformer_xrfi((float *)&db_out->block[curblock_out].data);

#if XRFI_TIME
		//gettimeofday(&tval_after2, NULL);
		//timersub(&tval_after2, &tval_before2, &tval_result2);
		//beam_time = (float)tval_result2.tv_usec/1000;
		clock_gettime(CLOCK_MONOTONIC, &tval_after2);
		time_taken2 = (float)(tval_after2.tv_sec - tval_before2.tv_sec)*1e6; // Time in seconds since epoch
		time_taken2 = time_taken2 + (float)(tval_after2.tv_nsec - tval_before2.tv_nsec)*1e-6; // Time in nanoseconds since 'tv_sec - start and end'
		beam_time = time_taken2;
#endif
                ///////////////////////////////////////////////////////////

		// Reinitialize projection matrix to identity matrix at the end of each iteration
		// This is similar to each time that a thread processes a block in HASHPIPE
        	proj_to_identity();

                //////////////////////////////////////////////////////////////////////////
                //>>>>    Timing correlator, projection matrix computation & RTBF    <<<<
                //////////////////////////////////////////////////////////////////////////
#if XRFI_TIME
                printf("run_correlator_xrfi()             -> time = %f ms\n", corr_time);
                printf("Computing projection matrix       -> time = %f ms\n", proj_time);
                printf("run_beamformer_xrfi()             -> time = %f ms\n", beam_time);
#endif
                //////////////////////////////////////////////////////////////////////////
           
                // Get block's starting mcnt for output block
                db_out->block[curblock_out].header.mcnt = db_in->block[curblock_in].header.mcnt_start;
                db_out->block[curblock_out].header.good_data = db_in->block[curblock_in].header.good_data;

#if VERBOSE
                printf("BF: Setting block %d, mcnt %lld as filled\n", curblock_out, (long long int)db_out->block[curblock_out].header.mcnt);
                printf("BF: Setting block %d, mcnt %lld as filled\n", curblock_out, (long long int)db_out->block[curblock_out].header.mcnt);
#endif
                // Mark output block as full and advance
                flag_gpu_beamformer_output_databuf_set_filled(db_out, curblock_out);
                curblock_out = (curblock_out + 1) % db_out->header.n_block;
                start_mcnt = last_mcnt + 1;
                last_mcnt = start_mcnt + Nm - 1;
		    
                // Mark input block as free
                flag_input_databuf_set_free(db_in, curblock_in);
                curblock_in = (curblock_in + 1) % db_in->header.n_block;
            }

        } else if (cur_state == CLEANUP) {
#if VERBOSE
            printf("RTBF: In Cleanup\n");
#endif

            hashpipe_status_lock_safe(&st);
            hgets(st.buf, "NETSTAT", 16, netstat);
            hashpipe_status_unlock_safe(&st);

            if (strcmp(netstat, "IDLE") == 0) {
                next_state = ACQUIRE;
                flag_databuf_clear((hashpipe_databuf_t *) db_out);
                printf("RTBF: Finished CLEANUP, clearing output databuf and returning to ACQUIRE\n");
            } else {
                next_state = CLEANUP;
                curblock_in = 0;
                curblock_out = 0;
                hashpipe_status_lock_safe(&st);
                hputl(st.buf, "CLEANA",1);
                hashpipe_status_unlock_safe(&st);
            }
        }

        // Next state processing
        hashpipe_status_lock_safe(&st);
        switch(next_state){
            case ACQUIRE: hputs(st.buf, status_key, "ACQUIRE"); break;
            case CLEANUP: hputs(st.buf, status_key, "CLEANUP"); break;        
        }
        hashpipe_status_unlock_safe(&st);
        cur_state = next_state;
        pthread_testcancel();
    }

    // Thread terminates after loop
    hashpipe_status_lock_busywait_safe(&st);
    // Unregister db_in->block[ii].data and db_out->block[ii].data
    // that were allocated in pinned memory with cudaHostUnregister()
    int ii;
    for (ii = 0; ii < N_INPUT_BLOCKS; ii++){
	unregister_data((signed char *)&db_in->block[ii].data);
        unregister_data((float *)&db_out->block[ii].data);
    } 
    unregister_data(c_data_out);
    printf("RTBF: Cleaning up gpu context...\n");
    corrCleanup();
    rtbfCleanup();
    hputs(st.buf, status_key, "terminated");
    hashpipe_status_unlock_safe(&st);
    return NULL;
}



// Thread description
static hashpipe_thread_desc_t b_thread = {
    name: "flag_xrfi_thread",
    skey: "BEAMSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {flag_input_databuf_create},
    obuf_desc: {flag_gpu_beamformer_output_databuf_create}
};

static __attribute__((constructor)) void ctor() {
    register_hashpipe_thread(&b_thread);
}

