
#ifndef _NQSORT_H_
#define _NQSORT_H_

#ifdef __cplusplus
extern "C" {
#endif

void nqsort(char* a, uint32_t n, uint32_t es, int32_t (*cmp)(char*, char*));

#ifdef __cplusplus
}
#endif

#endif
