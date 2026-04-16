#ifndef DAC_H
#define DAC_H

#include <stdint.h>

typedef struct {
    uint8_t **level;        // chunk arrays
    uint64_t **cont;        // continuation bitmaps
    uint32_t *lvl_size;     // size per level
    uint32_t nlevels;       // number of levels
    uint32_t n;             // number of elements
} dac_t;

/* Build DAC from array of uint32_t */
dac_t *dac_build(const uint32_t *A, uint32_t n);

/* Access i-th element */
uint32_t dac_access(const dac_t *D, uint32_t idx);

/* Free memory */
void dac_free(dac_t* D);

uint64_t dac_bits(const dac_t *D);

#endif /* DAC_H */
