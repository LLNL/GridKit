#include <FileIO.hpp>
#include <MatPowerUtils.hpp>
#include <Testing.hpp>
#include <iostream>

using namespace GridKit;
using namespace GridKit::Testing;
using namespace GridKit::MatPowerUtils;

namespace {

using IntT = int;
using RealT = double;

static const std::string matpower_data{
    R"(

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0   14	0;
	2	0	0	3	0   15	0;
	2	0	0	3	0   30	0;
	2	0	0	3	0   40	0;
	2	0	0	3	0   10	0;
];

)"};

} // namespace

int main(int argc, char **argv) {
  int fail = 0;
  std::vector<GenCostRow<IntT, RealT>> gencost_answer{
      {2, 0, 0, 3, {0, 14, 0}}, {2, 0, 0, 3, {0, 15, 0}},
      {2, 0, 0, 3, {0, 30, 0}}, {2, 0, 0, 3, {0, 40, 0}},
      {2, 0, 0, 3, {0, 10, 0}},
  };
  MatPower<IntT, RealT> mp;

  {
    std::istringstream iss(matpower_data);
    GridKit::readMatPower(mp, iss);
    if (!isEqual(mp.gencost, gencost_answer))
      fail++;
    std::cout << "After reading the gencost component, fail == " << fail
              << "\n";
  }

  std::cout << "Tests " << (fail ? "FAILED" : "PASSED") << "\n";

  return fail;
}
