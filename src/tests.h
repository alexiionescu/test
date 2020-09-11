#pragma once
void *test_genetics(void *user_data, const char *line, size_t size,FILE* out);
void ParseParams(char *input, int n, ...);
int ParseAllParams(char* input, int argc, char** argv);