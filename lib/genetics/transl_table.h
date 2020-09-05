#pragma once

extern char TRANSL_TABLE[64];
extern char STARTS_TABLE[64];
extern const char *TRANSL_TABLE_LONG[64];
extern char DNA_STRINGS[64][4];
extern char RNA_STRINGS[64][4];

void parse_transl_table(int n);