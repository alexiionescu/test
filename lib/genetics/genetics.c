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

#define CODON_1stBase 4
#define CODON_2ndBase 2
#define CODON_3rdBase 0
#define CODON_1stBaseMask 0x30
#define CODON_1stBaseMaskINV 0xCF
#define CODON_2ndBaseMask 0x0C
#define CODON_2ndBaseMaskINV 0xF3
#define CODON_3rdBaseMask 0x03
#define CODON_3rdBaseMaskINV 0xFC

#define GET_CODON_1stBase(codon) (codon >> CODON_1stBase & 0x3)
#define SET_CODON_1stBase(codon, bp) ((bp << CODON_1stBase) | (codon & CODON_1stBaseMaskINV))
#define GET_CODON_2ndBase(codon) (codon >> CODON_2ndBase & 0x3)
#define SET_CODON_2ndBase(codon, bp) ((bp << CODON_2ndBase) | (codon & CODON_2ndBaseMaskINV))
#define GET_CODON_3rdBase(codon) (codon & 0x03)
#define SET_CODON_3rdBase(codon, bp) (bp | (codon & CODON_3rdBaseMaskINV))

/**
 * @brief CODON ENCODING on a  byte 00<3rdBase><2ndBase><1stBase>
 * (U)T->0 00 
 *    C->1 01
 *    A->2 10
 *    G->3 11
 */
typedef char CODON_BASE;
// invert direction, exchange 1st and 3rd base
#define REVERSE_CODON(codon) ((codon & CODON_2ndBaseMask) | (codon >> CODON_1stBase & 0x3) | ((codon & 0x3) << CODON_1stBase))
// Get codon complement (U)T<->A, C<->G ( XOR 101010 )
#define COMPLEMENT_CODON(codon) (codon ^ 0x2A)

#define DNA_BUFFER_START 2
struct _GeneticsObj
{
    uint8_t *dna;
    DNA_DIR dnaDir;
    size_t dnaSize;
    uint8_t *dnaAllocBuffer;
    size_t dnaAllocSize;
    bool dnaInput;
    CODON_BASE lastCodonBase;
    FILE *out;
    uint8_t start_codon;
};

static int transl_table = 0;
/**
 * @brief Create a new genetics object.
 *        Use Genetics_Delete() to delete.
 * 
 * @return GeneticsObj* newly created object
 */
GeneticsObj *Genetics_New()
{
    if (0 == transl_table)
    {
        transl_table = 1;
        parse_transl_table(transl_table);
    }

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
    if (_this->start_codon == n)
        return;
    if (n > _this->start_codon)
    {
        for (int i = 0; i < _this->dnaSize; i++)
        {
            if (n - _this->start_codon == 1) //1->2,2->3
            {
                _this->dna[i - 1] <<= CODON_2ndBase;
                _this->dna[i - 1] &= 0x3f;
                _this->dna[i - 1] = SET_CODON_3rdBase(_this->dna[i - 1], GET_CODON_1stBase(_this->dna[i]));
            }
            else // 1->3
            {
                _this->dna[i - 1] <<= CODON_1stBase;
                _this->dna[i - 1] &= 0x3f;
                _this->dna[i - 1] = SET_CODON_2ndBase(_this->dna[i - 1], GET_CODON_1stBase(_this->dna[i]));
                _this->dna[i - 1] = SET_CODON_3rdBase(_this->dna[i - 1], GET_CODON_2ndBase(_this->dna[i]));
            }
        }

        if (n - _this->start_codon == 1) //1->2,2->3
        {
            _this->dna[_this->dnaSize - 1] <<= CODON_2ndBase;
            _this->dna[_this->dnaSize - 1] &= 0x3f;
        }
        else
        {
            _this->dna[_this->dnaSize - 1] <<= CODON_1stBase;
            _this->dna[_this->dnaSize - 1] &= 0x3f;
        }
        //  1  XX 01  10  11  XX 01' 10' 11'
        //  2  XX 10  11  01' XX 10' 11' 01"
        //  3  XX 11  01' 10' XX 11' 01" 10"
    }
    else
    {
        for (int r = _this->dnaSize - 1; r >= 0; r--)
        {
            if (_this->start_codon - n == 1) //3->2,2->1
            {
                _this->dna[r] >>= CODON_2ndBase;
                _this->dna[r] = SET_CODON_1stBase(_this->dna[r], GET_CODON_3rdBase(_this->dna[r - 1]));
            }
            else // 3->1
            {
                _this->dna[r] >>= CODON_1stBase;
                _this->dna[r] = SET_CODON_1stBase(_this->dna[r], GET_CODON_2ndBase(_this->dna[r - 1]));
                _this->dna[r] = SET_CODON_2ndBase(_this->dna[r], GET_CODON_3rdBase(_this->dna[r - 1]));
            }
        }

        if (_this->start_codon - n == 1) //3->2,2->1
        {
            _this->dna[-1] >>= CODON_2ndBase;
        }
        else // 3->1
        {
            _this->dna[-1] >>= CODON_1stBase;
        }
    }
    _this->start_codon = n;
}

/**
 * @brief Start DNA Input
 * 
 * @param _this genetics object
 * @param dir direction: 
 *              DNA_DIR_5_TO_3 : 5' to 3' strand
 *              DNA_DIR_3_TO_5 : 3' to 5' strand
 * @param code DNA code: a string containing letters A T(U) G C
 */
void Genetics_StartDNA(GeneticsObj *_this, DNA_DIR dir, const char *code)
{
    _this->dnaInput = true;
    if (_this->dnaAllocSize == 0)
    {
        _this->dnaAllocSize = 10240;
        _this->dnaAllocBuffer = (uint8_t *)malloc(_this->dnaAllocSize);
        _this->dna = _this->dnaAllocBuffer + DNA_BUFFER_START;
    }
    _this->dnaDir = dir;
    _this->dnaSize = 0;
    _this->dna[0] = 0;
    _this->lastCodonBase = CODON_1stBase;
    _this->start_codon = 1;
    Genetics_AddDNA(_this, code);
}

/**
 * @brief Add DNA Input
 * 
 * @param _this genetics object
 * @param code DNA code: a string containing letters A T(U) G C
 */
void Genetics_AddDNA(GeneticsObj *_this, const char *code)
{
    if (*code == ';' || *code == '>')
        return; //FASTA lines
    if (_this->dnaInput)
    {
        size_t codeSize = strlen(code);
        if (_this->dnaAllocSize <= _this->dnaSize + codeSize)
        {
            _this->dnaAllocSize = 2 * _this->dnaAllocSize + codeSize + 1;
            _this->dnaAllocBuffer = (uint8_t *)realloc(_this->dnaAllocBuffer, _this->dnaAllocSize);
            _this->dna = _this->dnaAllocBuffer + DNA_BUFFER_START;
        }
        while (*code)
        {
            switch (*code)
            {
            case 'T':
            case 't':
                _this->lastCodonBase -= 2;
                break;
            case 'U':
            case 'u':
                _this->lastCodonBase -= 2;
                break;
            case 'C':
            case 'c':
                _this->dna[_this->dnaSize] |= 1 << _this->lastCodonBase;
                _this->lastCodonBase -= 2;
                break;
            case 'A':
            case 'a':
                _this->dna[_this->dnaSize] |= 2 << _this->lastCodonBase;
                _this->lastCodonBase -= 2;
                break;
            case 'G':
            case 'g':
                _this->dna[_this->dnaSize] |= 3 << _this->lastCodonBase;
                _this->lastCodonBase -= 2;
                break;
            }
            if (_this->lastCodonBase < 0)
            {
                _this->lastCodonBase = CODON_1stBase;
                _this->dnaSize++;
                _this->dna[_this->dnaSize] = 0;
            }
            code++;
        }
    }
    else
    {
        fprintf(_this->out, "warning Genetics_AddDNA without DNA Start");
    }
}

/**
 * @brief Stop DNA Input
 * 
 * @param _this genetics object
 */
void Genetics_StopDNA(GeneticsObj *_this)
{
    if (_this->dnaInput)
    {
        _this->dnaInput = false;
        if (_this->lastCodonBase != CODON_1stBase)
        {
            _this->dnaSize++;
        }
    }
    else
    {
        fprintf(_this->out, "warning Genetics_StopDNA without DNA Start");
    }
}

/**
 * @brief Check DNA Input Status
 * 
 * @param _this genetics object
 * @return true if DNA Input is Started, false otherwise
 * 
 */
bool Genetics_DNAInput(GeneticsObj *_this)
{
    return _this->dnaInput;
}

#define PSTATE_NA -1
#define START_LINE_FMT  "\n%8lu "
static void PrintCodon(GeneticsObj *_this, size_t i, uint8_t codon,
                       DNA_PRINT_FlAGS flags, int *pstate)
{
    if (flags & DNA_PRINT_REVERSE)
        codon = REVERSE_CODON(codon);
    if (flags & DNA_PRINT_COMPLEMENT)
        codon = COMPLEMENT_CODON(codon);
    bool translChanged = false;
    if (flags & (DNA_PRINT_TRANSLATE | DNA_PRINT_TRANSLATE_LONG))
    {
        if (*pstate == PSTATE_NA && STARTS_TABLE[codon] == 'M')
        {
            if (flags & DNA_PRINT_TRANSLATE)
                fprintf(_this->out, START_LINE_FMT "M", i * 3 + _this->start_codon);
            else
                fprintf(_this->out, START_LINE_FMT "Met", i * 3 + _this->start_codon);
            translChanged = true;
            *pstate = i;
        }
        else if (*pstate != PSTATE_NA)
        {
            if (STARTS_TABLE[codon] == '*')
            {
                translChanged = true;
                *pstate = PSTATE_NA;
            }
            else if (STARTS_TABLE[codon] == 'M')
            {
                //fputc('>',_this->out); //alternative starts
            }
        }
    }

    if (flags & DNA_PRINT_TRANSLATE)
    {
        if (!translChanged && *pstate != PSTATE_NA && (i - *pstate) % 60 == 0)
            fprintf(_this->out, " ..." START_LINE_FMT, i * 3 + _this->start_codon);
    }
    else if (flags & DNA_PRINT_TRANSLATE_LONG)
    {
        if (!translChanged && *pstate != PSTATE_NA)
        {
            if ((i - *pstate) % 20 == 0)
                fprintf(_this->out, " ..." START_LINE_FMT, i * 3 + _this->start_codon);
            else
                fputc('-', _this->out);
        }
    }
    else
    {
        size_t offset = i;
        if (flags & DNA_PRINT_REVERSE) 
            offset = _this->dnaSize - i - 1;
        if (offset % 20 == 0)
        {
            size_t poffset = i * 3 + _this->start_codon;   
            if (flags & DNA_PRINT_REVERSE) 
                poffset += 2;
            fprintf(_this->out, START_LINE_FMT, poffset);
        }
        else
        {
            if (offset % 4 == 0)
                fputc(' ', _this->out);
        }
    }

    if (flags & DNA_PRINT_TRANSLATE)
    {
        if (!translChanged && *pstate != PSTATE_NA)
            fputc(TRANSL_TABLE[codon], _this->out);
    }
    else if (flags & DNA_PRINT_TRANSLATE_LONG)
    {
        if (!translChanged && *pstate != PSTATE_NA)
            fputs(TRANSL_TABLE_LONG[codon], _this->out);
    }
    else
    {
        if (flags & DNA_PRINT_RNA)
            fputs(RNA_STRINGS[codon], _this->out);
        else
            fputs(DNA_STRINGS[codon], _this->out);
    }
}

static void PrintLastCodon(GeneticsObj *_this, uint8_t codon, DNA_PRINT_FlAGS flags)
{
    if (flags & (DNA_PRINT_TRANSLATE | DNA_PRINT_TRANSLATE_LONG))
        return;

    int side = 2;
    int mid = 1;
    if (flags & DNA_PRINT_REVERSE)
    {
        codon = REVERSE_CODON(codon);
        side = 0;
    }
    if (flags & DNA_PRINT_COMPLEMENT)
        codon = COMPLEMENT_CODON(codon);

    char codonStr[4] = "aaa";
    if (flags & DNA_PRINT_RNA)
        memcpy(codonStr, RNA_STRINGS[codon], 3);
    else
        memcpy(codonStr, DNA_STRINGS[codon], 3);

    switch (_this->start_codon)
    {
    case 1:
        switch (_this->lastCodonBase)
        {
        case CODON_2ndBase:
            codonStr[mid] = '-';
        case CODON_3rdBase:
            codonStr[side] = '-';
        }
        break;
    case 2:
        switch (_this->lastCodonBase)
        {
        case CODON_3rdBase:
            codonStr[mid] = '-';
        case CODON_1stBase:
            codonStr[side] = '-';
        }
        break;
    case 3:
        switch (_this->lastCodonBase)
        {
        case CODON_1stBase:
            codonStr[mid] = '-';
        case CODON_2ndBase:
            codonStr[side] = '-';
        }
        break;
    }
    fputs(codonStr, _this->out);
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
        fprintf(_this->out, " /codon_start=%d", _this->start_codon);
    }
}

/**
 * @brief Print DNA Info
 * 
 * @param _this genetics object
 * @param flags \n 
 *          DNA_PRINT_REVERSE: Reverse order \n 
 *          DNA_PRINT_COMPLEMENT: Print Complement pairs \n 
 *          DNA_PRINT_RNA: DNA to RNA: T -> U \n 
 *          DNA_PRINT_TRANSLATE: Translate to protein (single letter)
 *          DNA_PRINT_TRANSLATE_LONG: Translate to protein (3 letters)
 */
void Genetics_PrintDNA(GeneticsObj *_this, DNA_PRINT_FlAGS flags)
{
    int pstate = PSTATE_NA;
    size_t printSize = _this->dnaSize;

    if (flags & DNA_PRINT_REVERSE)
    {
        PrintHeader(_this, true, flags);
        if (_this->start_codon != 1 || _this->lastCodonBase != CODON_1stBase)
        {
            printSize--;
            if (_this->start_codon == 3 && _this->lastCodonBase == CODON_2ndBase)
                printSize--;
            if (3 - _this->start_codon != _this->lastCodonBase / 2){
                fprintf(_this->out,START_LINE_FMT,(size_t)(printSize + 1) * 3 - _this->start_codon + 1);
                PrintLastCodon(_this, _this->dna[printSize], flags);
            }
        }
        for (int r = printSize - 1; r >= 0; r--)
        {
            PrintCodon(_this, r, _this->dna[r], flags, &pstate);
        }
        PrintHeader(_this, false, flags);
    }
    else
    {
        if (_this->start_codon != 1 || _this->lastCodonBase != CODON_1stBase)
        {
            printSize--;
            if (_this->start_codon == 3 && _this->lastCodonBase == CODON_2ndBase)
                printSize--;
        }
        PrintHeader(_this, true, flags);
        for (size_t i = 0; i < printSize; i++)
        {
            PrintCodon(_this, i, _this->dna[i], flags, &pstate);
        }
        if (_this->start_codon != 1 || _this->lastCodonBase != CODON_1stBase)
        {
            if (3 - _this->start_codon != _this->lastCodonBase / 2)
                PrintLastCodon(_this, _this->dna[printSize], flags);
        }
        PrintHeader(_this, false, flags);
    }
    fputs("\n-------------------------\n\n", _this->out);
}

/**
 * @brief print to current ouput file
 * 
 * @param _this genetics object
 * @param s string to print
 */
void Genetics_Print(GeneticsObj *_this, const char *s)
{
    fputs(s, _this->out);
}
/**
 * @brief Load a FASTA file
 * 
 * @param _this genetics object
 * @param filename   FASTA file name
 * @param search    string to search for in ^> lines. when found will load dna from next line to the next ^> line
 */
void Genetics_LoadFASTA(GeneticsObj *_this, const char *filename, const char *search)
{
    size_t len = 0;
    char *input = NULL;
    ssize_t insize;
    FILE *fInput = NULL;
    fInput = fopen(filename, "r");
    if (!fInput)
    {
        fprintf(stderr, "Error fopening fasta file '%s' : %s\n", filename, strerror(errno));
        return;
    }
    fprintf(_this->out, "Load FASTA file '%s' searching for '%s'\n", filename, search);
    bool bFound = false;
    Genetics_StartDNA(_this, DNA_DIR_5_TO_3, "");
    while (-1 != (insize = getline(&input, &len, fInput)))
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
        if (bFound) //dna line
            Genetics_AddDNA(_this, input);
    }
    Genetics_StopDNA(_this);
    fprintf(_this->out, "FASTA loaded. Found %lu bp.\n", _this->dnaSize * 3);
    fclose(fInput);
    free(input);
}