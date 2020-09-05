#include <config.h>

#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include <getopt.h>
#include <ctype.h>

#include "tests.h"

typedef void *(*TFunc)(void *user_data, const char *line);
struct test_menu
{
    const char *command;
    TFunc controler;
    void *user_data;
};

static struct test_menu test_menus[] = {
    {"genetics", test_genetics},
    {}};
static int menu_index = -1;

#define COLOR_OFF "\001\x1B[0m\002"
#define COLOR_RED "\001\x1B[0;91m\002"
#define COLOR_GREEN "\001\x1B[0;92m\002"
#define COLOR_YELLOW "\001\x1B[0;93m\002"
#define COLOR_BLUE "\001\x1B[0;94m\002"
#define COLOR_MAGENTA "\001\x1B[0;95m\002"
#define COLOR_BOLDGRAY "\001\x1B[1;30m\002"
#define COLOR_BOLDWHITE "\001\x1B[1;37m\002"
#define COLOR_HIGHLIGHT "\001\x1B[1;39m\002"

#define PROMPT_MAIN COLOR_BLUE "test" COLOR_OFF "$"
#define PROMPT_CMD COLOR_MAGENTA "%s" COLOR_OFF "$"

bool ProcessNewInput(int line_no, bool fromFile, char *input, size_t insize)
{
    size_t term = insize - 1;
    while (input[term] == '\r' || input[term] == '\n')
    {
        input[term] = 0;
        term--;
    }
    char *line = input;
    while (isspace(*line))
        line++;

    if (*line == '#')
        return true;

    if (!strncasecmp(line, "quit", 4))
    {
        if (!fromFile)
            puts("Bye!");
        return false;
    }

    if (menu_index != -1)
    {
        if (!strncasecmp(line, "back", 4) || !strncasecmp(line, "exit", 4))
        {
            menu_index = -1;
            if (!fromFile)
                fputs(PROMPT_MAIN, stdout);
            return true;
        }
    }

    if (menu_index == -1)
    {
        if (!strncasecmp(line, "echo", 4))
        {
            line += 4;
            while (isspace(*line))
                line++;
            puts(line);
            if (!fromFile)
                fputs(PROMPT_MAIN, stdout);
            return true;
        }

        for (int i = 0; test_menus[i].command; i++)
        {
            if (!strcasecmp(line, test_menus[i].command))
            {
                menu_index = i;
                break;
            }
        }
    }

    if (menu_index == -1)
    {
        if (fromFile)
            printf("Unknown menu:%d: %s\n", line_no, line);
        else
            puts("Unknown menu" PROMPT_MAIN);
        return true;
    }

    if (fromFile)
    {
        test_menus[menu_index].user_data = test_menus[menu_index].controler(test_menus[menu_index].user_data, line);
    }
    else
    {
        test_menus[menu_index].user_data = test_menus[menu_index].controler(test_menus[menu_index].user_data, line);
        printf(PROMPT_CMD, test_menus[menu_index].command);
    }

    return true;
}

int main(int argc, char **argv)
{

    const char *UsagePrint = "Usage: "
                             "[ -f | --filename INPUTFILE ]\n";
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"filename", required_argument, 0, 'f'},
        {}};
    int opt;
    FILE *fInput = NULL;
    while (-1 != (opt = getopt_long(argc, argv, "hf:", long_options, NULL)))
    {
        switch (opt)
        {
        case 'h':
            puts(UsagePrint);
            exit(0);
            break;
        case 'f':
        {
            char *filename = optarg;
            while (isspace(*filename))
                filename++;
            fInput = fopen(filename, "r");
            if (!fInput)
                fprintf(stderr, "Error fopening input file '%s' : %s\n", optarg, strerror(errno));
        }
        break;
        case '?':
            break;
        }
    }

    size_t len = 0;
    char *input = NULL;
    ssize_t insize;
    int line_no = 1;
    bool quit = false;
    if (fInput)
    {
        while (-1 != (insize = getline(&input, &len, fInput)))
        {
            if (!ProcessNewInput(line_no++, true, input, insize))
            {
                quit = true;
                break;
            }
        }
        fclose(fInput);
    }
    if (!quit)
    {
        if (menu_index == -1)
            fputs(PROMPT_MAIN, stdout);
        else
            printf(PROMPT_CMD, test_menus[menu_index].command);

        while (-1 != (insize = getline(&input, &len, stdin)))
        {
            if (!ProcessNewInput(line_no++, false, input, insize))
            {
                quit = true;
                break;
            }
        }
    }

    free(input);
    for (int i = 0; test_menus[i].command; i++)
    {
        if (test_menus[menu_index].user_data)
            test_menus[menu_index].controler(test_menus[menu_index].user_data, NULL);
    }
    return 0;
}
