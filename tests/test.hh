#ifndef TEST_H
#define TEST_H

#include <iostream>

class Test
{
private:
    void run();
public:
    int tests = 10000;
    void end_test(int err)
    {
        if (err)
            std::cout << err << " error";
        else
            std::cout << "CLEAR";
        std::cout << std::endl;
    }
};

#endif
