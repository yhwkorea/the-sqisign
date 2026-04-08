#ifndef BITSIZE_TRACKER_H
#define BITSIZE_TRACKER_H

/**
 * @file bitsize_tracker.h
 * @brief Track maximum intermediate bit size during computation
 *
 * Used for benchmarking HNF vs MLLL intermediate value growth.
 * Call tracker_reset() before, tracker_get_max() after.
 * Inside computation, call tracker_update() with each intermediate ibz_t.
 */

#include <quaternion.h>

#ifdef BITSIZE_TRACKER_ENABLE
/* Global tracker state (not thread-safe, but fine for benchmarks).
 * Defined as weak so that any TU can link without a separate definition.
 * mlll_benchmark.c provides strong definitions when present. */
#ifdef __GNUC__
__attribute__((weak)) int _bitsize_tracker_max = 0;
__attribute__((weak)) int _bitsize_tracker_enabled = 0;
#else
extern int _bitsize_tracker_max;
extern int _bitsize_tracker_enabled;
#endif

static inline void tracker_reset(void) {
    _bitsize_tracker_max = 0;
    _bitsize_tracker_enabled = 1;
}

static inline void tracker_disable(void) {
    _bitsize_tracker_enabled = 0;
}

static inline int tracker_get_max(void) {
    return _bitsize_tracker_max;
}

static inline void tracker_update_ibz(const ibz_t *v) {
    if (_bitsize_tracker_enabled) {
        int b = ibz_bitsize(v);
        if (b > _bitsize_tracker_max)
            _bitsize_tracker_max = b;
    }
}

static inline void tracker_update_vec4(const ibz_vec_4_t *v) {
    if (_bitsize_tracker_enabled) {
        for (int i = 0; i < 4; i++)
            tracker_update_ibz(&((*v)[i]));
    }
}

static inline void tracker_update_mat4x4(const ibz_mat_4x4_t *m) {
    if (_bitsize_tracker_enabled) {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                tracker_update_ibz(&((*m)[i][j]));
    }
}
#else
/* No-op stubs when tracker is disabled */
static inline void tracker_reset(void) {}
static inline void tracker_disable(void) {}
static inline int tracker_get_max(void) { return 0; }
static inline void tracker_update_ibz(const ibz_t *v) { (void)v; }
static inline void tracker_update_vec4(const ibz_vec_4_t *v) { (void)v; }
static inline void tracker_update_mat4x4(const ibz_mat_4x4_t *m) { (void)m; }
#endif

#endif
