#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <signal.h>

static long fail_at = -1;  /* fail at this allocation count */
static long call_count = 0;

static void init_once(void) {
    const char *s;
    static int inited = 0;
    if (inited) return;
    inited = 1;
    s = getenv("FAIL_MALLOC_AT");
    if (s) {
        char *end = NULL;
        long v = strtol(s, &end, 10);
        if (end != s && v > 0) fail_at = v;
    }
}

void * __real_malloc(size_t);
void * __real_calloc(size_t, size_t);


static int should_fail(void){
    init_once();
    call_count++;
    return (fail_at > 0 && call_count == fail_at);
}

void * __wrap_malloc(size_t sz){
    if (should_fail()){
        fprintf(stderr, "[wrap] malloc failing at call #%ld\n", call_count);
        errno = ENOMEM;
        return NULL;
    }
    return __real_malloc(sz);
}

void * __wrap_calloc(size_t n, size_t sz) {
    if (should_fail()) {
        fprintf(stderr, "[wrap] calloc failing at call #%ld\n", call_count);
        errno = ENOMEM;
        return NULL;
    }
    return __real_calloc(n, sz);
}