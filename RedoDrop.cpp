//
//  main.cpp
//  DropletShape
//
//  Created by Kara Maki on 9/16/21.
//  Copyright © 2021 Kara Maki. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string>
#include "Redo.h"

int main(int argc, const char* argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";

    Drop Try;
    Try.Initialization();
    std::cout << "Initilization done\n";
    Try.CalculateTimeArrays();
    Try.printCurrentMassInformation("Xitwocompleted41");
    return 0;
}
