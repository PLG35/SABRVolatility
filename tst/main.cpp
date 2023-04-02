#include "test_smile_curve.h"
#include "test_planner_v1.h"
#include "test_planner_v2.h"

using namespace std;

int main(int argc, char *argv[])
{

    MainTestPlan myTests = MainTestPlan();
    myTests.run();

/*
 * TEST OF COPY
 *  MainTestPlan myTests2 = myTests;
 *  myTests2.run();
 */

    // launch_tests();

    return 0;
}
