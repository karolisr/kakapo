#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "fqgzlib.h"


int main(int argc, char const *argv[])
{
    if(argc == 2) {
        int avg_read_len;
        avg_read_len = gz_fq_avg_read_len((char *)argv[1]);
        printf("%s: %d\n", "Average read length", avg_read_len);
    }
    return 0;
}
