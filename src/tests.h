#pragma once

void *test_genetics(void *user_data, const char *line, size_t size,FILE* out);
void ParseParams(char *input, int n, ...);
int ParseAllParams(char* input, int argc, char** argv);

#define COLOR_OFF "\001\x1B[0m\002"
#define COLOR_RED "\001\x1B[0;91m\002"
#define COLOR_GREEN "\001\x1B[0;92m\002"
#define COLOR_YELLOW "\001\x1B[0;93m\002"
#define COLOR_BLUE "\001\x1B[0;94m\002"
#define COLOR_MAGENTA "\001\x1B[0;95m\002"
#define COLOR_BOLDGRAY "\001\x1B[1;30m\002"
#define COLOR_BOLDWHITE "\001\x1B[1;37m\002"
#define COLOR_HIGHLIGHT "\001\x1B[1;39m\002"
#define HELP_START_LINE "\n\t\t"
typedef struct _menu_item {
    const char* command;
    const char* params;
    const char* description;
}menu_help_item;
extern menu_help_item MenuShared[];
extern menu_help_item MenuGlobal[];
void PrintHelp(menu_help_item* menu);
bool PrintMenuHelp(const char* line, menu_help_item* menu_help);
