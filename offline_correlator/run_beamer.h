#ifndef __RUN_BEAMER_H
#define __RUN_BEAMER_H

#define TIMESTEPS_PER_CALL 100
// ?
int run_beamer(int8_t* device_data,
    int8_t *data,
    uint8_t * device_results,
    uint8_t *results,
    size_t step_data_size,
    size_t step_results_size
);

#endif
