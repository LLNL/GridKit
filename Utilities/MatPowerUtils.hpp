#pragma once
#include <regex>
#include <string>
#include <iomanip>
#include <sstream>

/**
 *
 * @file Utilities/MatPowerUtils.hpp
 * @author Asher Mancinelli <asher.mancinelli@pnnl.gov>
 *
 * @remark `std::stringstream` is preferred over `operator+(std::string, ...)`
 * due to stringstream's lack of reallocation on append.
 *
 */

namespace GridKit
{

namespace MatPowerUtils
{

  template <typename IntT = int, typename RealT = double>
  struct BusRow {
    IntT bus_i;
    IntT type;
    IntT Pd;
    RealT Qd;
    IntT Gs;
    IntT Bs;
    IntT area;
    RealT Vm;
    RealT Va;
    IntT baseKV;
    IntT zone;
    RealT Vmax;
    RealT Vmin;

    inline std::string str() const
    {
      std::stringstream ss;
      std::cerr << std::setw(10) << bus_i << std::setw(10) << type << std::setw(10)
         << Pd << std::setw(10) << Qd << std::setw(10) << Gs << std::setw(10)
         << Bs << std::setw(10) << area << std::setw(10) << Vm << std::setw(10)
         << Va << std::setw(10) << baseKV << std::setw(10) << zone
         << std::setw(10) << Vmax << std::setw(10) << Vmin;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename IntT = int, typename RealT = double>
  struct GenRow {
    IntT bus;
    RealT Pg;
    RealT Qg;
    RealT Qmax;
    RealT Qmin;
    RealT Vg;
    IntT mBase;
    IntT status;
    IntT Pmax;
    IntT Pmin;
    IntT Pc1;
    IntT Pc2;
    IntT Qc1min;
    IntT Qc1max;
    IntT Qc2min;
    IntT Qc2max;
    IntT ramp_agc;
    IntT ramp_10;
    IntT ramp_30;
    IntT ramp_q;
    IntT apf;

    inline std::string str() const
    {
      std::stringstream ss;
      ss << std::setw(10) << bus << std::setw(10) << Pg << std::setw(10) << Qg
         << std::setw(10) << Qmax << std::setw(10) << Qmin << std::setw(10)
         << Vg << std::setw(10) << mBase << std::setw(10) << status
         << std::setw(10) << Pmax << std::setw(10) << Pmin << std::setw(10)
         << Pc1 << std::setw(10) << Pc2 << std::setw(10) << Qc1min
         << std::setw(10) << Qc1max << std::setw(10) << Qc2min << std::setw(10)
         << Qc2max << std::setw(10) << ramp_agc << std::setw(10) << ramp_10
         << std::setw(10) << ramp_30 << std::setw(10) << ramp_q << std::setw(10)
         << apf;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename IntT = int, typename RealT = double>
  struct BranchRow {
    IntT fbus;
    IntT tbus;
    RealT r;
    RealT x;
    RealT b;
    IntT rateA;
    IntT rateB;
    IntT rateC;
    IntT ratio;
    IntT angle;
    IntT status;
    IntT angmin;
    IntT angmax;

    inline std::string str() const
    {
      std::stringstream ss;
      ss << std::setw(10) << fbus << std::setw(10) << tbus << std::setw(10) << r
         << std::setw(10) << x << std::setw(10) << b << std::setw(10) << rateA
         << std::setw(10) << rateB << std::setw(10) << rateC << std::setw(10)
         << ratio << std::setw(10) << angle << std::setw(10) << status
         << std::setw(10) << angmin << std::setw(10) << angmax;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename IntT = int, typename RealT = double>
  struct GenCostRow {
    IntT kind;
    IntT startup;
    IntT shutdown;
    IntT n;
    std::vector<RealT> rest;

    inline std::string str() const
    {
      std::stringstream ss;
      ss << std::setw(10) << kind << std::setw(10) << startup << std::setw(10)
         << shutdown << std::setw(10) << n;
      for (const auto& val : rest)
        ss << std::setw(10) << val;
      ss << "\n";
      return ss.str();
    }
  };

  template <typename IntT = int, typename RealT = double>
  struct MatPower {
    using BusRowT = BusRow<IntT, RealT>;
    using GenRowT = GenRow<IntT, RealT>;
    using BranchRowT = BranchRow<IntT, RealT>;
    using GenCostRowT = GenCostRow<IntT, RealT>;

    std::string version;
    IntT baseMVA;
    std::vector<BusRowT> bus;
    std::vector<GenRowT> gen;
    std::vector<BranchRowT> branch;
    std::vector<GenCostRowT> gencost;

    // Not sure if these should be in this struct... Not all matpower files
    // I found contained them.
    //
    // std::string name;
    // std::vector<std::string> bus_name;

    inline std::string str() const
    {
      std::stringstream ss;
      ss << "Version: " << version << "\n"
         << "Base MVA: " << baseMVA << "\n";

      ss << "Bus:\n";
      for (const auto& v : bus)
        ss << bus.str();

      ss << "Gen:\n";
      for (const auto& v : gen)
        ss << gen.str();

      ss << "Branch:\n";
      for (const auto& v : branch)
        ss << branch.str();

      ss << "GenCost:\n";
      for (const auto& v : gencost)
        ss << gencost.str();

      ss << "\n";

      return ss.str();
    }
  };

}  // namespace MatPower
}  // namespace GridKit
