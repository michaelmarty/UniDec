#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "isogenpep.h"


int main(const int argc, char *argv[])
{
    printf("Testing IsoGen from exe\n");
    char* filename = argv[1];
    printf("Testing IsoGen on %s\n", filename);

    run_file(filename);
}
