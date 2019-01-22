//
//  main.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include <iostream>
#include <unistd.h>
#include "test.h"
#include "utility/timer.h"
#include "utility/logger.h"

using namespace std;

int main(int argc, const char* argv[])
{
    LOGGER_CONF("feyncalc.log", "main", Logger::file_on | Logger::screen_on, INFO, INFO);
    RunTest();
    return 0;
}
