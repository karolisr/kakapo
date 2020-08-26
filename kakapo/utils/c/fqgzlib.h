ssize_t gzgetdelim(char **buf, size_t *bufsiz, int delimiter, gzFile fp);
ssize_t gzgetline(char **buf, size_t *bufsiz, gzFile fp);
unsigned int gz_fq_avg_read_len(char *path);
