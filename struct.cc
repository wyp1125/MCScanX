/* Helper functions for error msgs for allocating memory, some are from
 * UCSC Jim Kent's library
 */

#include "struct.h"

void progress(const char *format, ...)
/* Print progress message */
{
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    fprintf(stdout, "\n");
    va_end(args);
}

void err(const char *format, ...)
/* Print error message but do not exit */
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Error] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

void warn(const char *format, ...)
/* Print error message but do not exit */
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Warn] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
}

void errAbort(const char *format, ...)
/* Print error message and exit */
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[Error] ");
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
    exit(1);
}

long clock1000()
/* A millisecond clock. */
{
    struct timeval tv;
    static long origSec;
    gettimeofday(&tv, NULL);
    if (origSec == 0) origSec = tv.tv_sec;
    return (tv.tv_sec-origSec)*1000 + tv.tv_usec / 1000;
}

void uglyTime(const char *label, ...)
/* Print label and how long it's been since last call.  Call with
 * a NULL label to initialize. */
{
    static long lastTime = 0;
    long time = clock1000();
    va_list args;
    va_start(args, label);
    if (label != NULL)
    {
        vfprintf(stdout, label, args);
        fprintf(stdout, " [%.3f seconds elapsed]\n", (time - lastTime)/1000.);
    }
    lastTime = time;
    va_end(args);
}

FILE *mustOpen(const char *fileName, const char *mode)
/* Open a file or die */
{
    FILE *f;
    /* gcc 4.2 has "deprecated conversion from string constant" problem"*/
    char *modeName = (char *)"";

    if (sameString(fileName, "stdin")) return stdin;
    if (sameString(fileName, "stdout")) return stdout;
    if ((f = fopen(fileName, mode)) == NULL)
    {
        if (mode)
        {
            if (mode[0] == 'r') modeName = (char *)" to read";
            else if (mode[0] == 'w') modeName = (char *)" to write";
            else if (mode[0] == 'a') modeName = (char *)" to append";
        }
        errAbort("Can't open %s%s: %s", fileName, modeName, strerror(errno));
    }
    return f;
}

