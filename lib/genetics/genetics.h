#pragma once
typedef struct _GeneticsObj GeneticsObj;
GeneticsObj *Genetics_New();
void Genetics_Delete(GeneticsObj *_this);

#define DNA_DIR_5_TO_3 0
#define DNA_DIR_3_TO_5 1
typedef char DNA_DIR;

size_t Genetics_StartDNA(GeneticsObj *_this, DNA_DIR dir, const char *code);
size_t Genetics_AddDNA(GeneticsObj *_this, const char *code);
void Genetics_StopDNA(GeneticsObj *_this);
bool Genetics_DNAInput(GeneticsObj *_this);
void Genetics_LoadFASTA(GeneticsObj *_this, size_t start,size_t stop, const char *filename, const char *search);

#define DNA_PRINT_REVERSE 0x0001
#define DNA_PRINT_COMPLEMENT 0x0002
#define DNA_PRINT_RNA 0x0004
#define DNA_PRINT_TRANSLATE 0x0008
#define DNA_PRINT_TRANSLATE_LONG 0x0010
#define DNA_PRINT_TRANSLATE_CORRELATE 0x0020
typedef uint32_t DNA_PRINT_FlAGS;

void Genetics_Print(GeneticsObj *_this, const char *s);
void Genetics_PrintDNA(GeneticsObj *_this, DNA_PRINT_FlAGS flags);
void Genetics_SetOutput(GeneticsObj *_this, FILE *out);
void Genetics_SetCodonStart(GeneticsObj *_this, int n);