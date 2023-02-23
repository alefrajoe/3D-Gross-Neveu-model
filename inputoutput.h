#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H

#include "macro.h"
#include "info.h"

// input output functions
bool ReadCorrectInput(char **argv);
bool CreateAllRequiredDirectories(void);
void StringObservables(char * string_observables, const struct Info const * info);
void StringAcceptance(char * string_acceptance, const struct Info const * info);
void StringTerma(char * string_terma, char * string_terma_copy, const struct Info const * info);
void OpenAndWriteObservablesTitle(FILE * opf, char * string_title);
void OpenAndWriteAcceptanceTitle(FILE * opf, char * string_title, const struct Info const * info);

// analysis function
double Mean(const double const * array, int num_elem);
#endif