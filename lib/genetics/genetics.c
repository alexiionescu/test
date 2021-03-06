#include <config.h>

#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>

#include "genetics.h"
#include "transl_table.h"


/**
 * @brief Base Pair Encoding (1 byte per pair)
 * (U)T->0 00 
 *    C->1 01
 *    A->2 10
 *    G->3 11
 */
#define CODON(b1,b2,b3) (((b1 & 0x3)<<4)|((b2 & 0x3)<<2)|(b3 & 0x3))
#define COMPLEMENT_CODON(codon) (codon ^ 0x2A)  // XOR 101010   T <-> A, C <-> G
#define COMPLEMENT(b) (b ^ 0x2)                 // XOR     10   T <-> A, C <-> G

#define DNA_BUFFER_START 0
struct _GeneticsObj
{
    uint8_t *dna;
    DNA_DIR dnaDir;
    size_t dnaSize;
    uint8_t *dnaAllocBuffer;
    size_t dnaAllocSize;
    bool dnaInput;
    FILE *out;
    uint8_t start_codon;
    size_t inputFileOffset;
    bool fileBegin;
    size_t* spliceData;
    int spliceSize;
};

static uint8_t START_CODON[3];
/**
 * @brief static function that setup translation table
 * 
 * @param n translation table > 0
 */
void Genetics_SetTranslationTable(int n)
{
    static int transl_table = 0;
    if (n > 0 && n != transl_table)
    {
        transl_table = n;
        int s = parse_transl_table(transl_table);
        if(s == 0)
        {
            fputs("Translation tables have no start codons!!!\n",stderr);
        }
        else
        {
            START_CODON[0] = (s & (0x3 << 4)) >> 4;
            START_CODON[1] = (s & (0x3 << 2)) >> 2;
            START_CODON[2] = s & 0x3;
        }
    }
}


/**
 * @brief Create a new genetics object.
 *        Use Genetics_Delete() to delete.
 * 
 * @return GeneticsObj* newly created object
 */
GeneticsObj *Genetics_New()
{
    Genetics_SetTranslationTable(1);
    GeneticsObj *_this = (GeneticsObj *)malloc(sizeof(GeneticsObj));
    memset(_this, 0, sizeof(GeneticsObj));
    _this->out = stdout;
    return _this;
}

/**
 * @brief Delete a genetics object
 * 
 * @param _this pointer to a genetics object created with Genetics_New()
 */
void Genetics_Delete(GeneticsObj *_this)
{
    if (_this->dnaAllocBuffer)
        free(_this->dnaAllocBuffer);
    free(_this);
}

/**
 * @brief set output file
 * 
 * @param _this pointer to a genetics object created with Genetics_New()
 * @param out filestream for output
 */
void Genetics_SetOutput(GeneticsObj *_this, FILE *out)
{
    _this->out = out;
}

/**
 * @brief Set start read frame base pair
 * 
 * @param _this 
 * @param n start codon 1(default), 2 or 3
 */
void Genetics_SetCodonStart(GeneticsObj *_this, int n)
{
    if(n >= 1 && n < _this->dnaSize)
        _this->start_codon = n;
}

/**
 * @brief Find Start Codon
 * 
 * @param _this genetics object
 * @param flags \n 
 *          DNA_PRINT_REVERSE: Find on Reverse strand; \n 
 */
bool Genetics_FindStart(GeneticsObj *_this, DNA_PRINT_FlAGS flags)
{
    uint8_t s1,s2,s3;
    if(flags&DNA_PRINT_COMPLEMENT)
    {
        s1 = COMPLEMENT(START_CODON[0]);   
        s2 = COMPLEMENT(START_CODON[1]);   
        s3 = COMPLEMENT(START_CODON[2]); 
    }
    else
    {
        s1 = START_CODON[0];   
        s2 = START_CODON[1];   
        s3 = START_CODON[2];   
    }
     
    if (_this->dnaDir == DNA_DIR_3_TO_5)
    {
        flags ^= DNA_PRINT_REVERSE;
    }

    if (flags & DNA_PRINT_REVERSE)
    {
        size_t poffset = _this->inputFileOffset + _this->dnaSize;
        int splice = _this->spliceSize - 1;
        while(_this->spliceData && splice > 0 && poffset < _this->spliceData[splice]) 
            splice--;
        for (int r = _this->dnaSize; r >= 2; r--, poffset --)
        {
            uint8_t b1, b2, b3;
            b1 = _this->dna[r];
            b2 = _this->dna[r - 1];
            b3 = _this->dna[r - 2];
            int cut = 0;
            if (_this->spliceData && splice > 0)
            {
                int s = _this->spliceData[splice] - poffset;
                if (s >= -1) // could be -1, 0, 1
                {
                    cut = _this->spliceData[splice] - _this->spliceData[splice - 1];
                    if (s == 1)
                    {
                        if (r - cut - 2 < 0)
                            break;
                        b1 = _this->dna[r - cut];
                        b2 = _this->dna[r - cut - 1];
                        b3 = _this->dna[r - cut - 2];
                        if (s1 == b1 && s2 == b2 && s3 == b3)
                        {
                            _this->start_codon = _this->dnaSize - r - 2;
                            return true;
                        }
                        r -= cut - 2;
                        poffset -= cut - 2;
                        splice -= 2;
                        continue;
                    }
                    if (s == 0)
                    {
                        if (r - cut < 0)
                            break;
                        b2 = _this->dna[r - cut];
                    }
                    if (r - cut - 1 < 0)
                        break;
                    b3 = _this->dna[r - cut - 1];
                }
            }
            if (s1 == b1 && s2 == b2 && s3 == b3)
            {
                _this->start_codon = _this->dnaSize - r;
                return true;
            }
        }
    }
    else
    {
        size_t poffset = _this->inputFileOffset + 1;
        int splice = 0;
        while(_this->spliceData && splice < _this->spliceSize && poffset > _this->spliceData[splice]) 
            splice++;
        for (int i = 0; i < _this->dnaSize - 2; i++, poffset++)
        {
            uint8_t b1, b2, b3;
            b1 = _this->dna[i];
            b2 = _this->dna[i + 1];
            b3 = _this->dna[i + 2];
            int cut = 0;
            if (_this->spliceData && splice < _this->spliceSize)
            {
                int s = poffset - _this->spliceData[splice];
                if (s >= -1) // could be -1, 0, 1
                {
                    cut = _this->spliceData[splice + 1] - _this->spliceData[splice];
                    if (s == 1)
                    {
                        if (i + cut + 1 >= _this->dnaSize)
                            break;
                        b1 = _this->dna[i + cut -1];
                        b2 = _this->dna[i + cut];
                        b3 = _this->dna[i + cut + 1];
                        if (s1 == b1 && s2 == b2 && s3 == b3)
                        {
                            _this->start_codon = i-2;
                            return true;
                        }
                        i += cut - 2;
                        poffset += cut - 2;
                        splice += 2;
                        continue;
                    }
                    if (s == 0)
                    {
                        if (i + cut >= _this->dnaSize)
                            break;
                        b2 = _this->dna[i + cut];
                    }
                    if (i + 1 + cut >= _this->dnaSize)
                        break;
                    b3 = _this->dna[i + 1 + cut];
                }
            }
            if (s1 == b1 && s2 == b2 && s3 == b3)
            {
                _this->start_codon = i + 1;
                return true;
            }
        }
    }

    return false;
}

/**
 * @brief Start DNA Input
 * 
 * @param _this genetics object
 * @param dir direction: 
 *              DNA_DIR_5_TO_3 : 5' to 3' strand
 *              DNA_DIR_3_TO_5 : 3' to 5' strand
 * @param code DNA code: a string containing letters A T(U) G C
 * 
 * @return number of bp added
 */
size_t Genetics_StartDNA(GeneticsObj *_this, DNA_DIR dir, const char *code)
{
    _this->dnaInput = true;
    if (_this->dnaAllocSize == 0)
    {
        _this->dnaAllocSize = 102400;
        _this->dnaAllocBuffer = (uint8_t *)malloc(_this->dnaAllocSize);
        _this->dna = _this->dnaAllocBuffer + DNA_BUFFER_START;
    }
    _this->dnaDir = dir;
    _this->dnaSize = 0;
    _this->dna[0] = 0;
    _this->start_codon = 1;
    _this->inputFileOffset = 0;
    _this->fileBegin = true;
    return Genetics_AddDNA(_this, code);
}

/**
 * @brief Add DNA Input
 * 
 * @param _this genetics object
 * @param code DNA code: a string containing letters A T(U) G C
 * 
 * @return number of bp added
 */
size_t Genetics_AddDNA(GeneticsObj *_this, const char *code)
{
    if (*code == ';' || *code == '>')
        return 0; //FASTA lines
    size_t bp = 0;
    if (_this->dnaInput)
    {
        size_t codeSize = strlen(code);
        if (_this->dnaAllocSize <= _this->dnaSize + codeSize)
        {
            _this->dnaAllocSize = 10 * _this->dnaAllocSize;
            _this->dnaAllocBuffer = (uint8_t *)realloc(_this->dnaAllocBuffer, _this->dnaAllocSize);
            _this->dna = _this->dnaAllocBuffer + DNA_BUFFER_START;
        }
        while (*code)
        {
            switch (*code)
            {
            case 'T':
            case 't':
                _this->dna[_this->dnaSize++] = 0;
                if(_this->fileBegin) _this->fileBegin = false; 
                bp++;
                break;
            case 'U':
            case 'u':
                _this->dna[_this->dnaSize++] = 0;
                if(_this->fileBegin) _this->fileBegin = false; 
                bp++;
                break;
            case 'C':
            case 'c':
                _this->dna[_this->dnaSize++] = 1;
                if(_this->fileBegin) _this->fileBegin = false; 
                bp++;
                break;
            case 'A':
            case 'a':
                _this->dna[_this->dnaSize++] = 2;
                if(_this->fileBegin) _this->fileBegin = false; 
                bp++;
                break;
            case 'G':
            case 'g':
                _this->dna[_this->dnaSize++] = 3;
                if(_this->fileBegin) _this->fileBegin = false; 
                bp++;
                break;
            case 'N':
            case 'n':
                if(_this->fileBegin) _this->inputFileOffset++;
                break;
            }
            code++;
        }
    }
    else
    {
        fprintf(_this->out, "warning Genetics_AddDNA without DNA Start");
    }
    return bp;
}

/**
 * @brief Stop DNA Input
 * 
 * @param _this genetics object
 */
void Genetics_StopDNA(GeneticsObj *_this)
{
    _this->dnaInput = false;
}

/**
 * @brief Check DNA Input Status
 * 
 * @param _this genetics object
 * @return DNA_DIR_NONE if not started, direction DNA_DIR_5_TO_3 or DNA_DIR_3_TO_5 otherwise
 * 
 */
int Genetics_DNAInput(GeneticsObj *_this)
{
    return _this->dnaInput ? _this->dnaDir : DNA_DIR_NONE;
}

#define PSTATE_NA -1
#define START_LINE_FMT      "\n%9lu "
#define START_LINE_EMPTY "\n          " //same spaces as START_LINE_FMT empty 
#define CODONS_PER_LINE 20
#define PROTEINS_BP_PER_LINE 210
#define PROTEINS_LONG_BP_PER_LINE 60
static void PrintCodon(FILE* out, size_t bufferOffset, uint8_t b1, uint8_t b2, uint8_t b3,
                       DNA_PRINT_FlAGS flags, int *pstate, size_t poffset, size_t printOffset)
{
    uint8_t codon = CODON(b1,b2,b3);
    if (flags & DNA_PRINT_COMPLEMENT)
        codon = COMPLEMENT_CODON(codon);
    bool translChanged = false;

    if (flags & (DNA_PRINT_TRANSLATE | DNA_PRINT_TRANSLATE_LONG))
    {
        if (*pstate == PSTATE_NA && STARTS_TABLE[codon] == 'M')
        {
            if(flags&DNA_PRINT_TRANSLATE_CORRELATE)
            {
                if (flags & DNA_PRINT_TRANSLATE)
                    fputs("M   ",out);
                else
                    fputs("Met-",out);
            }
            else
            {
                if (flags & DNA_PRINT_TRANSLATE)
                    fprintf(out, START_LINE_FMT "M", poffset);
                else
                    fprintf(out, START_LINE_FMT "Met-", poffset);    
            }
            translChanged = true;
            *pstate = bufferOffset;
        }
        else if (*pstate != PSTATE_NA)
        {
            if (STARTS_TABLE[codon] == '*')
            {
                translChanged = true;
                *pstate = PSTATE_NA;
            }
        }
    }

    if (flags & DNA_PRINT_TRANSLATE)
    {
        if(!(flags&DNA_PRINT_TRANSLATE_CORRELATE))
        {
            if (!translChanged && *pstate != PSTATE_NA && printOffset % PROTEINS_BP_PER_LINE == 0)
                fprintf(out, " ..." START_LINE_FMT, poffset);
        }
    }
    else if (flags & DNA_PRINT_TRANSLATE_LONG)
    {
        if(!(flags&DNA_PRINT_TRANSLATE_CORRELATE))
        {
            if (!translChanged && *pstate != PSTATE_NA && printOffset % PROTEINS_LONG_BP_PER_LINE == 0)
                fprintf(out, " ..." START_LINE_FMT, poffset);                
        }        
    }
    else
    {
        if (printOffset % CODONS_PER_LINE == 0)
        {
            fprintf(out, START_LINE_FMT, poffset);
        }
        else 
        {
            fputc(' ', out);
        }
    }

    if (flags & DNA_PRINT_TRANSLATE)
    {
        if (!translChanged){ 
            if(*pstate != PSTATE_NA)
            {
                fputc(TRANSL_TABLE[codon], out);
                if(flags&DNA_PRINT_TRANSLATE_CORRELATE)
                    fputs("   ",out);
            }
            else if(flags&DNA_PRINT_TRANSLATE_CORRELATE)
                fputs("    ",out);
        }
        else if(translChanged && *pstate == PSTATE_NA && flags&DNA_PRINT_TRANSLATE_CORRELATE)
            fputs("    ",out);
    }
    else if (flags & DNA_PRINT_TRANSLATE_LONG)
    {
        if (!translChanged)
        { 
            if(*pstate != PSTATE_NA)
            {
                fputs(TRANSL_TABLE_LONG[codon], out);
                fputc('-',out);
            }
            else if(flags&DNA_PRINT_TRANSLATE_CORRELATE)
                fputs("    ",out);
        }
        else if(translChanged && *pstate == PSTATE_NA && flags&DNA_PRINT_TRANSLATE_CORRELATE)
            fputs("    ",out);
    }
    else
    {
        if (flags & DNA_PRINT_RNA)
            fputs(RNA_STRINGS[codon], out);
        else
            fputs(DNA_STRINGS[codon], out);
    }
}

static void PrintHeader(GeneticsObj *_this, bool begin, DNA_PRINT_FlAGS flags)
{
    if (flags & (DNA_PRINT_TRANSLATE | DNA_PRINT_TRANSLATE_LONG))
    {
        if (begin)
            fputs("\nNH2", _this->out);
        else
            fputs("\nCOOH", _this->out);
    }
    else
    {
        if (((flags & DNA_PRINT_COMPLEMENT) && (flags & DNA_PRINT_REVERSE)) ||
            (!(flags & DNA_PRINT_COMPLEMENT) && !(flags & DNA_PRINT_REVERSE)))
        { //same order
            if ((begin && _this->dnaDir == DNA_DIR_5_TO_3) || (!begin && _this->dnaDir == DNA_DIR_3_TO_5))
                fputs("\n5'", _this->out);
            else
                fputs("\n3'", _this->out);
        }
        else
        { // reverse order
            if ((begin && _this->dnaDir == DNA_DIR_5_TO_3) || (!begin && _this->dnaDir == DNA_DIR_3_TO_5))
                fputs("\n3'", _this->out);
            else
                fputs("\n5'", _this->out);
        }
    }
    if (begin && _this->start_codon != 1)
    {
        if(flags & DNA_PRINT_REVERSE)
            fprintf(_this->out, " /codon_start <%d (%lu)", _this->start_codon, _this->inputFileOffset + 1 +_this->dnaSize - _this->start_codon);
        else
            fprintf(_this->out, " /codon_start %d> (%lu)", _this->start_codon, _this->inputFileOffset + _this->start_codon);
    }
}

#define END_PRINT_STRING    "\n-------------------------\n\n"
/**
 * @brief Print DNA Info
 * 
 * @param _this genetics object
 * @param flags \n 
 *          DNA_PRINT_REVERSE: Reverse order; \n 
 *          DNA_PRINT_COMPLEMENT: Print Complement pairs; \n 
 *          DNA_PRINT_RNA: DNA to RNA; \n 
 *          DNA_PRINT_TRANSLATE: Translate to protein (single letter); \n
 *          DNA_PRINT_TRANSLATE_LONG: Translate to protein (3 letters); \n
 *          DNA_PRINT_TRANSLATE_CORRELATE: Show BP and Translate;
 */
void Genetics_PrintDNA(GeneticsObj *_this, DNA_PRINT_FlAGS flags)
{
    int pstate = PSTATE_NA;
    if(_this->dnaSize == 0) {
        PrintHeader(_this, true, flags);
        PrintHeader(_this, false, flags);
        fputs(END_PRINT_STRING, _this->out);
        return;
    }
    if(_this->dnaDir == DNA_DIR_3_TO_5)
    {
        flags ^= DNA_PRINT_REVERSE;
    }

    if(flags&DNA_PRINT_TRANSLATE_CORRELATE && _this->spliceData){
        fprintf(stderr,"ERROR: Splice is not supported with correlate translation\n");
        return;
    }

    size_t printOffset = 0;
    
    bool printCorrelation = false;
    bool endCorrelation = false;
    DNA_PRINT_FlAGS pflags = flags;
    if(flags&DNA_PRINT_TRANSLATE_CORRELATE) 
        pflags = flags & ~(DNA_PRINT_TRANSLATE_LONG|DNA_PRINT_TRANSLATE);
    if (flags & DNA_PRINT_REVERSE)
    {
        PrintHeader(_this, true, flags);
        size_t poffset = _this->inputFileOffset + 1 +_this->dnaSize - _this->start_codon;   
        size_t cor_r = _this->dnaSize - _this->start_codon, cor_poffset = poffset;
        int splice = _this->spliceSize - 1;
        while(_this->spliceData && splice > 0 && poffset < _this->spliceData[splice]) 
            splice--;
        for (int r = cor_r; r >= 2; r-=3, poffset -= 3, printOffset++)
        {
            if(flags&DNA_PRINT_TRANSLATE_CORRELATE && printOffset > 0 && printOffset % CODONS_PER_LINE == 0 && !endCorrelation) 
            {
                printCorrelation = !printCorrelation;
                if(printCorrelation)
                {
                    r = cor_r;
                    poffset = cor_poffset; 
                    fputs(START_LINE_EMPTY,_this->out);
                    pflags = flags;
                }
                else
                {
                    cor_r = r;
                    cor_poffset = poffset;
                    pflags = flags & ~(DNA_PRINT_TRANSLATE_LONG|DNA_PRINT_TRANSLATE);
                }
            }
            uint8_t b1,b2,b3;
            b1 = _this->dna[r];
            b2 = _this->dna[r-1];
            b3 = _this->dna[r-2];
            int cut = 0;
            if(_this->spliceData && splice > 0)
            {
                int s = _this->spliceData[splice] - poffset;
                if(s >= -1) // could be -1, 0, 1
                {
                    cut = _this->spliceData[splice] - _this->spliceData[splice-1];
                    splice -= 2;
                    if(s==1) {
                        r -= cut - 4; 
                        poffset -= cut - 4;
                        printOffset = -1;
                        continue;
                    }
                    if(s==0)
                    {
                        if(r - cut < 0) break;
                        b2 = _this->dna[r - cut]; 
                    }
                    r--;
                    poffset--; 
                    if(r - cut < 0) break; 
                    b3 = _this->dna[r - cut];  
                    r--;
                    poffset--;
                }                
            }

            PrintCodon(_this->out, r, b1,b2,b3, pflags, &pstate, poffset, printOffset);

            if(cut > 0)
            {
                r -= cut - 3; 
                poffset -= cut - 3;
                printOffset = -1;
            }

            if(flags&DNA_PRINT_TRANSLATE_CORRELATE && !printCorrelation && r - 3 < 2) 
            {
                printCorrelation = true;
                endCorrelation = true;
                r = cor_r + 3;
                poffset = cor_poffset + 3; 
                fputs(START_LINE_EMPTY,_this->out);
                pflags = flags;
            }
        }
        PrintHeader(_this, false, flags);
    }
    else
    {
        PrintHeader(_this, true, flags);
        size_t poffset = _this->inputFileOffset + _this->start_codon;
        size_t cor_i = _this->start_codon - 1, cor_poffset = poffset;
        int splice = 0;
        while(_this->spliceData && splice < _this->spliceSize && poffset > _this->spliceData[splice]) 
            splice++;
        for (size_t i = cor_i; i < _this->dnaSize - 2; i+=3, poffset += 3, printOffset++)
        {
            if(flags&DNA_PRINT_TRANSLATE_CORRELATE && printOffset > 0 && printOffset % CODONS_PER_LINE == 0 && !endCorrelation) 
            {
                printCorrelation = !printCorrelation;
                if(printCorrelation)
                {
                    i = cor_i;
                    poffset = cor_poffset; 
                    fputs(START_LINE_EMPTY,_this->out);
                    pflags = flags;
                }
                else
                {
                    cor_i = i;
                    cor_poffset = poffset;
                    pflags = flags & ~(DNA_PRINT_TRANSLATE_LONG|DNA_PRINT_TRANSLATE);
                }
            }
            uint8_t b1,b2,b3;
            b1 = _this->dna[i];
            b2 = _this->dna[i+1];
            b3 = _this->dna[i+2];
            int cut = 0;
            if(_this->spliceData && splice < _this->spliceSize)
            {
                int s = poffset - _this->spliceData[splice];
                if(s >= -1) // could be -1, 0, 1
                {
                    cut = _this->spliceData[splice+1] - _this->spliceData[splice];
                    splice += 2;
                    if(s==1) {
                        i += cut - 4; 
                        poffset += cut - 4;
                        printOffset = -1;
                        continue;
                    }
                    if(s==0)
                    {
                        if(i + cut >= _this->dnaSize) break;
                        b2 = _this->dna[i + cut];   
                    }
                    i++;
                    poffset++;
                    if(i + cut >= _this->dnaSize) break;
                    b3 = _this->dna[i + cut];
                    i++;
                    poffset++;
                }                
            }
            
            PrintCodon(_this->out, i, b1,b2,b3, pflags, &pstate, poffset, printOffset);
            
            if(cut > 0)
            {
                i += cut - 3; 
                poffset += cut - 3;
                printOffset = -1;
            }

            if(flags&DNA_PRINT_TRANSLATE_CORRELATE && !printCorrelation && i + 3 >= _this->dnaSize - 2) 
            {
                printCorrelation = true;
                endCorrelation = true;
                fputs(START_LINE_EMPTY,_this->out);
                i = cor_i - 3;
                poffset = cor_poffset - 3; 
                pflags = flags;
            }
        }
        PrintHeader(_this, false, flags);
    }
    fputs(END_PRINT_STRING, _this->out);
}

/**
 * @brief Load a FASTA file
 * 
 * @param _this genetics object
 * @param start dna bp start
 * @param stop  dna bp stop
 * @param filename   FASTA file name
 * @param search    string to search for in ^> lines. when found will load dna from next line to the next ^> line
 */
void Genetics_LoadFASTA(GeneticsObj *_this, size_t start, size_t stop, const char *filename, const char *search)
{
    size_t len = 0;
    char *input = NULL;
    ssize_t insize;
    FILE *fInput = NULL;
    fInput = fopen(filename, "r");
    size_t lines = 0;
    if (!fInput)
    {
        fprintf(stderr, "Error fopening fasta file '%s' : %s\n", filename, strerror(errno));
        return;
    }
    if(start > 0 && stop <= start)
    {
        fprintf(stderr, "Error fopening fasta file '%s' : stop %lu is less then start %lu\n", filename, stop,start);
        return;
    }
    fprintf(_this->out, "Load FASTA file '%s' searching for '%s'\n", filename, search);
    bool bFound = false;
    Genetics_StartDNA(_this, DNA_DIR_5_TO_3, "");
    size_t read;
    while (-1 != (read = (insize = getline(&input, &len, fInput))))
    {
        if (*input == ';')
        { //fasta commented line
            continue;
        }
        if (*input == '>')
        { //fasta genome description line
            bFound = *search == 0 || NULL != strstr(input, search);
            if (bFound)
                fputs(input, _this->out);
            continue;
        }
        if (bFound)
        {
            if(start > 0){
                if(_this->inputFileOffset == 0)
                {
                    start--;
                    _this->inputFileOffset = start;
                    fseek(fInput, _this->inputFileOffset + _this->inputFileOffset / (read - 1) - read,SEEK_CUR);
                    continue;
                }
                else if(start + read - 1 > stop)
                {
                    input[stop - start] = 0;
                }   
            }
            start += Genetics_AddDNA(_this, input);
            lines++;

            if(start == stop)
            {
                break;
            }
        }
    }
    Genetics_StopDNA(_this);
    fprintf(_this->out, "FASTA loaded. Found %lu bp on %lu lines.\n", _this->dnaSize, lines);
    fclose(fInput);
    free(input);
}

/**
 * @brief Splice Data for next print
 * 
 * @param _this genetics object
 * @param n     size of data
 * @param data  data (splice offsets)
 */
void Genetics_Splice(GeneticsObj *_this, int n, size_t* data)
{
    if(n % 2 == 1)
    {
        fprintf(stderr,"ERROR: Splice data must have even size\n");
        return;
    }
    if(n == 0)
    {
        if(_this->spliceData) 
            free(_this->spliceData);
        _this->spliceSize = 0;
    }
    else
    {
        if(_this->spliceData) 
            free(_this->spliceData);
        _this->spliceSize = 0;
        _this->spliceData = malloc(n * sizeof(size_t));
        for(int i=0; i<n;i++)
        {
            int k=0;
            for(;k<_this->spliceSize; k++)
            {
                if(_this->spliceData[k] > data[i]) break;
            }
            if(_this->spliceSize > k)
                memmove(_this->spliceData + k + 1, _this->spliceData + k,(_this->spliceSize - k) * sizeof(size_t));
            _this->spliceData[k] = data[i];
            _this->spliceSize++;
        }
        
    }
}