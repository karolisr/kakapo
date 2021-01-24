#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "gzlib.h"


ssize_t gzgetdelim(char **buf, size_t *bufsiz, int delimiter, gzFile fp) {
    char *ptr, *eptr;


    if (*buf == NULL || *bufsiz == 0) {
        *bufsiz = BUFSIZ;
        if ((*buf = malloc(*bufsiz)) == NULL)
            return -1;
    }

    for (ptr = *buf, eptr = *buf + *bufsiz;;) {
        char c = gzgetc(fp);
        if (c == -1) {
            if (gzeof(fp)) {
                ssize_t diff = (ssize_t) (ptr - *buf);
                if (diff != 0) {
                    *ptr = '\0';
                    return diff;
                }
            }
            return -1;
        }
        *ptr++ = c;
        if (c == delimiter) {
            *ptr = '\0';
            return ptr - *buf;
        }
        if (ptr + 2 >= eptr) {
            char *nbuf;
            size_t nbufsiz = *bufsiz * 2;
            ssize_t d = ptr - *buf;
            if ((nbuf = realloc(*buf, nbufsiz)) == NULL)
                return -1;
            *buf = nbuf;
            *bufsiz = nbufsiz;
            eptr = nbuf + nbufsiz;
            ptr = nbuf + d;
        }
    }
}


ssize_t gzgetline(char **buf, size_t *bufsiz, gzFile fp) {
    return gzgetdelim(buf, bufsiz, '\n', fp);
}
