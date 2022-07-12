#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "gzlib.h"
#include "fqlib.h"


unsigned long long int fq_avg_read_len(char *path) {
    FILE *stream;
    char *line = NULL;
    size_t len = 0;
    unsigned int read;
    unsigned long long int totalSeqLen = 0;
    unsigned int totalReads;
    unsigned long int totalLines = 0;

    stream = fopen(path, "r");
    if (stream == NULL)
        exit(EXIT_FAILURE);

    while ((read = getline(&line, &len, stream)) != -1) {
        totalLines++;
        if (totalLines % 4 == 2) {
            if (read - 1 > 1023) {
                return 1023;
            }
            totalSeqLen += read - 1;
        }
    }

    free(line);
    fclose(stream);
    totalReads = totalLines / 4;
    return totalSeqLen / totalReads;
}


unsigned long long int fq_avg_read_len_gz(char *path) {
    gzFile stream;
    char *line = NULL;
    size_t len = 0;
    unsigned int read;
    unsigned long long int totalSeqLen = 0;
    unsigned int totalReads;
    unsigned long int totalLines = 0;

    stream = gzopen(path, "r");
    if (stream == NULL)
        exit(EXIT_FAILURE);

    if (gzbuffer(stream, 1024 * 1024) == -1) {
        exit(EXIT_FAILURE);
    }

    while ((read = gzgetline(&line, &len, stream)) != -1) {
        totalLines++;
        if (totalLines % 4 == 2) {
            if (read - 1 > 1023) {
                return 1023;
            }
            totalSeqLen += read - 1;
        }
    }

    free(line);
    gzclose(stream);
    totalReads = totalLines / 4;
    return totalSeqLen / totalReads;
}
