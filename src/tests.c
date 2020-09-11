#include <config.h>

#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <stdarg.h>

#include "tests.h"

#include "lib/genetics/genetics.h"

DNA_PRINT_FlAGS GetPrintFlag(char *sp, char *ep)
{
    if (!strncasecmp("translate", sp, ep - sp))
        return DNA_PRINT_TRANSLATE;
    if (!strncasecmp("translate_long", sp, ep - sp))
        return DNA_PRINT_TRANSLATE_LONG;
    if (!strncasecmp("compl", sp, ep - sp))
        return DNA_PRINT_COMPLEMENT;
    if (!strncasecmp("rev", sp, ep - sp))
        return DNA_PRINT_REVERSE;
    if (!strncasecmp("rna", sp, ep - sp))
        return DNA_PRINT_RNA;
    if (!strncasecmp("cor", sp, ep - sp))
        return DNA_PRINT_TRANSLATE_CORRELATE;  
    return 0;
}

void *test_genetics(void *user_data, const char *line, size_t size, FILE* out)
{
    if (!user_data)
    { //init
        return Genetics_New();
    }
    if (line == NULL)
    { //cleanup
        Genetics_Delete(user_data);
        return NULL;
    }
    Genetics_SetOutput(user_data, out);
    if (!strncasecmp("splice", line, 6))
    {
        char* params[100];
        static const int psize = sizeof(params)/sizeof(char*);
        size_t spliceData[psize];
        int n = ParseAllParams((char *)line + 6, psize, params);
        for(int i = 0;i<n; i++)
        {
            spliceData[i] = strtoul(params[i],NULL,10);
        }
        Genetics_Splice(user_data, n,spliceData);
        return user_data;
    }
    if (!strncasecmp("load_fasta", line, 10))
    {
        char *filename, *search,*start,*stop;
        ParseParams((char *)line + 10, 4, &start,&stop, &filename, &search);
        Genetics_LoadFASTA(user_data, strtoul(start,NULL,10), strtoul(stop,NULL,10), filename, search);
        return user_data;
    }
    if (!strncasecmp("print", line, 5))
    {
        DNA_PRINT_FlAGS flags = 0;
        char *sp = (char *)line + 5;
        char *ep = sp;
        while (*ep)
        {
            if (isspace(*ep))
            {
                if (sp != ep)
                {
                    flags |= GetPrintFlag(sp, ep);
                    sp = ep;
                }
                sp++;
            }
            ep++;
        }
        if (sp != ep)
            flags |= GetPrintFlag(sp, ep);
        Genetics_PrintDNA(user_data, flags);
        return user_data;
    }
    if (!strncmp("5'", line, 2))
    {
        if (Genetics_DNAInput(user_data))
        {
            Genetics_StopDNA(user_data);
        }
        else
        {
            Genetics_StartDNA(user_data, DNA_DIR_5_TO_3, line + 2);
            if (line[size - 1] == '\'' && line[size - 2] == '3')
                Genetics_StopDNA(user_data);
        }
        return user_data;
    }
    if (!strncmp("3'", line, 2))
    {
        if (Genetics_DNAInput(user_data))
        {
            Genetics_StopDNA(user_data);
        }
        else
        {
            Genetics_StartDNA(user_data, DNA_DIR_3_TO_5, line + 2);
            if (line[size - 1] == '\'' && line[size - 2] == '5')
                Genetics_StopDNA(user_data);
        }
        return user_data;
    }
    if (!strncmp("codon_start", line, 11))
    {
        char *codon_start;
        ParseParams((char *)line + 11, 1, &codon_start);
        Genetics_SetCodonStart(user_data, atoi(codon_start));
        return user_data;
    }

    if (Genetics_DNAInput(user_data))
    {
        Genetics_AddDNA(user_data, line);
    }
    return user_data;
}

void ParseParams(char *input, int n, ...)
{
    va_list args;
    va_start(args, n);
    for (int i = 0; i < n - 1; i++)
    {
        char **p = va_arg(args, char **);
        while (*input && isspace(*input))
            input++;
        *p = input;
        while (*input && !isspace(*input))
            input++;
        if (*input)
        {
            *input = 0;
            input++;
        }
    }
    char **p = va_arg(args, char **);
    while (*input && isspace(*input))
        input++;
    *p = input;
    va_end(args);
}

int ParseAllParams(char* input, int argc, char** argv)
{
    int n = 0;
    while(*input && n < argc)
    {
        while (*input && isspace(*input))
            input++;
        argv[n++] = input;  
        while (*input && !isspace(*input))
            input++;
        if (*input)
        {
            *input = 0;
            input++;
        }
    }
    return n;
}