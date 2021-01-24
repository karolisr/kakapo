#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "fqlib.h"


int main(int argc, char const *argv[])
{
    if(argc == 2) {
        int avg_read_len;
        avg_read_len = fq_avg_read_len_gz((char *)argv[1]);
        printf("%s: %d\n", "Average read length", avg_read_len);
    }
    return 0;
}
