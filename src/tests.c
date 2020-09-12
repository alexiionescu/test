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

menu_help_item MenuGenetics[] =
{
    { "5' | 3'", "[...]" , "start/end dna sequence."
            HELP_START_LINE "Everything between this 2 commands is considered dna sequence."
            HELP_START_LINE "Can be used on the same line for example 5'agtaaggcc3'."},
    { "load_fasta", "start stop filename [search]" , "load fasta file from start to stop offset."
            HELP_START_LINE "Use 0 for start/stop to load all."
            HELP_START_LINE "Option <search> option will search for fasta > lines and if found will start from next line."},
    { "splice", "[s1 s2 s3 s4 ...]" , "splice dna sequence based on exons boundaries"
            HELP_START_LINE "exons are [start s1] [s2 s3] ... [sN stop]"
            HELP_START_LINE "introns are [s1+1 s2-1] [s3+1 s4-1] ..."},
    { "print", "[print flags]", "print dna sequence. Flags are:"
            HELP_START_LINE "\t translate : translate to proteins (single letters)"
            HELP_START_LINE "\t translate_long : translate to proteins (3 letters)"
            HELP_START_LINE "\t compl : complement dna sequence"
            HELP_START_LINE "\t codon_start : set codon reading frame (default 1)"
            HELP_START_LINE "\t rev : reverse dna sequence. Use with 'compl' and translate for reverse strand translation."
            HELP_START_LINE "\t cor : use with translate to show dna sequence and translation correlated. Does not work with splice."
            HELP_START_LINE "\t rna : print rna instead of dna (T becomes U)"
            },
    {}
};

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
    if(PrintMenuHelp(line,MenuGenetics)) return user_data;

    int dir;
    if (DNA_DIR_NONE != (dir = Genetics_DNAInput(user_data)))
    {
        Genetics_AddDNA(user_data, line);
        if(dir == DNA_DIR_5_TO_3)
        {
            if (line[size - 1] == '\'' && line[size - 2] == '3')
                Genetics_StopDNA(user_data);
        }
        else
        {
            if (line[size - 1] == '\'' && line[size - 2] == '5')
                Genetics_StopDNA(user_data);
        }
        
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

void PrintHelp(menu_help_item* menu)
{
    while(menu->command)
    {
        printf(COLOR_HIGHLIGHT "%s" COLOR_OFF " %s" HELP_START_LINE "%s\n",menu->command, menu->params, menu->description);
        menu++;
    }
}

bool PrintMenuHelp(const char* line, menu_help_item* menu_help)
{
    if (!strncasecmp(line, "help", 4))
    {
        PrintHelp(menu_help);
        puts("-------");
        PrintHelp(MenuShared);
        PrintHelp(MenuGlobal);
        return true;
    }

    return false;
}