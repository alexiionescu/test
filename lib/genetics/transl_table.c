#include <config.h>

#include <string.h>
#include <ctype.h>

#include "transl_table.h"

// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
const char *transl_tables[] = {
    " AAs   = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    "Starts = ---m------**--*----m---------------M----------------------------"
    "Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    "Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    "Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"};

const char *AAS_LONG[] = {
    "Ala", //A Alanine
    "---", //B ----
    "Cys", //C Cysteine
    "Asp", //D Aspartic Acid
    "Glu", //E Glutamic Acid
    "Phe", //F Phenylalanine
    "Gly", //G Glycine
    "His", //H Histidine
    "Ile", //I Isoleucine
    "---", //J -----
    "Lys", //K Lysine
    "Leu", //L Leucine
    "Met", //M Methionine
    "Asn", //N Asparagine
    "---", //O -----
    "Pro", //P Proline
    "Gln", //Q Glutamine
    "Arg", //R Arginine
    "Ser", //S Serine
    "Thr", //T Threonine
    "---", //U -----
    "Val", //V Valine
    "Trp", //W Tryptophan
    "---", //X -----
    "Tyr", //Y Tyrosine
    "---"  //Z -----
};

char TRANSL_TABLE[64];
char STARTS_TABLE[64];
const char *TRANSL_TABLE_LONG[64];
char DNA_STRINGS[64][4];
char RNA_STRINGS[64][4];

void parse_transl_table(int n)
{
    const char *tbl = transl_tables[n - 1];
    char *aas = 2 + strchr(tbl, '=');
    char *starts = 2 + strchr(aas, '=');
    char *base1 = 2 + strchr(starts, '=');
    char *base2 = 2 + strchr(base1, '=');
    char *base3 = 2 + strchr(base2, '=');

    for (int i = 0; i < 64; i++)
    {
        STARTS_TABLE[i] = starts[i];
        TRANSL_TABLE[i] = aas[i];
        TRANSL_TABLE_LONG[i] = AAS_LONG[aas[i] - 'A'];
        if (base1[i] == 'T')
        {
            DNA_STRINGS[i][0] = 't';
            RNA_STRINGS[i][0] = 'u';
        }
        else
        {
            DNA_STRINGS[i][0] = tolower(base1[i]);
            RNA_STRINGS[i][0] = tolower(base1[i]);
        }

        if (base2[i] == 'T')
        {
            DNA_STRINGS[i][1] = 't';
            RNA_STRINGS[i][1] = 'u';
        }
        else
        {
            DNA_STRINGS[i][1] = tolower(base2[i]);
            RNA_STRINGS[i][1] = tolower(base2[i]);
        }

        if (base3[i] == 'T')
        {
            DNA_STRINGS[i][2] = 't';
            RNA_STRINGS[i][2] = 'u';
        }
        else
        {
            DNA_STRINGS[i][2] = tolower(base3[i]);
            RNA_STRINGS[i][2] = tolower(base3[i]);
        }
        DNA_STRINGS[i][3] = 0;
        RNA_STRINGS[i][3] = 0;
    }
}