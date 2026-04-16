#include "dac.h"
#include <stdlib.h>
#include <string.h>

#define CHUNK_BITS 7
#define CHUNK_MASK ((1U << CHUNK_BITS) - 1)

#define popc(x) __builtin_popcountll(x)

/* -------------------- RANK1 -------------------- */
static inline uint32_t rank1(const uint64_t *bm, uint32_t i) {
    uint32_t w = i >> 6;
    uint32_t off = i & 63;
    uint32_t r = 0;

    for (uint32_t k = 0; k < w; k++)
        r += popc(bm[k]);

    if (off)
        r += popc(bm[w] & ((1ULL << off) - 1));

    return r;
}

/* -------------------- BUILD -------------------- */
dac_t *dac_build(const uint32_t *A, uint32_t n) {
    dac_t *D = (dac_t*) malloc(sizeof(dac_t));
    D->n = n;

    uint32_t *cur = (uint32_t*) malloc(sizeof(uint32_t)*n);
    memcpy(cur, A, sizeof(uint32_t)*n);

    D->level = malloc(sizeof(uint8_t*) * 32);
    D->cont  = malloc(sizeof(uint64_t*) * 32);
    D->lvl_size = malloc(sizeof(uint32_t) * 32);

    uint32_t active = n;
    uint32_t lvl = 0;

    while (active) {
        uint8_t *lev = aligned_alloc(64, active);
        uint64_t *bm = aligned_alloc(64, ((active + 63) >> 6) * sizeof(uint64_t));
        memset(bm, 0, ((active + 63) >> 6) * sizeof(uint64_t));

        uint32_t next_active = 0;

        for (uint32_t i = 0; i < active; i++) {
            lev[i] = cur[i] & CHUNK_MASK;
            uint32_t more = cur[i] >> CHUNK_BITS;

            if (more) {
                bm[i >> 6] |= 1ULL << (i & 63);
                cur[next_active++] = more;
            }
        }

        D->level[lvl] = lev;
        D->cont[lvl]  = bm;
        D->lvl_size[lvl] = active;

        active = next_active;
        lvl++;
    }

    D->nlevels = lvl;
    free(cur);
    return D;
}

/* -------------------- ACCESS -------------------- */
uint32_t dac_access(const dac_t *D, uint32_t idx) {
    uint32_t value = 0;
    uint32_t shift = 0;

    for (uint32_t lvl = 0; lvl < D->nlevels; lvl++) {
        value |= (uint32_t)D->level[lvl][idx] << shift;

        uint64_t mask = 1ULL << (idx & 63);
        if (!(D->cont[lvl][idx >> 6] & mask))
            break;

        idx = rank1(D->cont[lvl], idx);
        shift += CHUNK_BITS;
    }

    return value;
}

/* -------------------- FREE -------------------- */
void dac_free(dac_t *D) {
    if (!D) return;

    for (uint32_t i = 0; i < D->nlevels; i++) {
        free(D->level[i]);
        free(D->cont[i]);
    }

    free(D->level);
    free(D->cont);
    free(D->lvl_size);
    free(D);
}

uint64_t dac_bits(const dac_t *D) {
    if (!D) return 0;

    uint64_t bits = 0;

    /* size of level arrays (each entry is 1 byte = 8 bits) */
    for (uint32_t lvl = 0; lvl < D->nlevels; lvl++) {
        bits += (uint64_t)D->lvl_size[lvl] * 8;
    }

    /* size of continuation bitmaps */
    for (uint32_t lvl = 0; lvl < D->nlevels; lvl++) {
        uint32_t words = (D->lvl_size[lvl] + 63) >> 6;  // number of 64-bit words
        bits += (uint64_t)words * 64;
    }

    /* pointer arrays overhead (optional but accurate) */
    bits += (uint64_t)D->nlevels * (
        sizeof(uint8_t*) * 8 +
        sizeof(uint64_t*) * 8 +
        sizeof(uint32_t) * 8
    );

    /* structure metadata */
    bits += 8 * (
        sizeof(dac_t)          // struct itself
    );

    return bits;
}
