#include "test_planner_v2.h"
#include "test_planner_calibrator.h"

using namespace std;

int main(int argc, char *argv[])
{
    // PlanUnitaryTestSmile myTestsSmile = PlanUnitaryTestSmile();
    // myTestsSmile.run();

    double forward = 0.01;
    int calDate = 0;
    int optMat = 365;
    vector<double> strikes = {0.01, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05};
    vector<double> volatilities = {150, 120, 105, 100, 103, 110, 120};
    TestCalibrator myTestsCalibrator = TestCalibrator();
    myTestsCalibrator.fillStructure(&strikes, calDate, optMat);
    myTestsCalibrator.fillMarketData(forward, &volatilities);
    myTestsCalibrator.run();

    return 0;
}
