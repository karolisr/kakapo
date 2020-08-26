#include <stdio.h>
#include <stdlib.h>

#include "fqlib.h"

unsigned int fq_avg_read_len(char *path) {
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
            totalSeqLen += read - 1;
        }
    }

    free(line);
    fclose(stream);
    totalReads = totalLines / 4;
    return totalSeqLen / totalReads;
}
