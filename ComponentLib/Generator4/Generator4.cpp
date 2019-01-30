/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit. For details, see github.com/LLNL/GridKit
 * Please also read the LICENSE file.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the disclaimer below.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the disclaimer (as noted below) in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the LLNS/LLNL nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
 * SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Lawrence Livermore National Laboratory is operated by Lawrence Livermore
 * National Security, LLC, for the U.S. Department of Energy, National
 * Nuclear Security Administration under Contract DE-AC52-07NA27344.
 *
 * This document was prepared as an account of work sponsored by an agency
 * of the United States government. Neither the United States government nor
 * Lawrence Livermore National Security, LLC, nor any of their employees
 * makes any warranty, expressed or implied, or assumes any legal liability
 * or responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents that
 * its use would not infringe privately owned rights. Reference herein to
 * any specific commercial product, process, or service by trade name,
 * trademark, manufacturer, or otherwise does not necessarily constitute or
 * imply its endorsement, recommendation, or favoring by the United States
 * government or Lawrence Livermore National Security, LLC. The views and
 * opinions of authors expressed herein do not necessarily state or reflect
 * those of the United States government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */


#include <iostream>
#include <cmath>
#include <ComponentLib/Bus/BaseBus.hpp>
#include "Generator4.hpp"

namespace ModelLib {

/*!
 * @brief Constructor for a simple generator model
 *
 * Arguments passed to ModelEvaluatorImpl:
 * - Number of equations = 4 differential + 2 algebraic = 6
 * - Number of quadratures = 1
 * - Number of optimization parameters = 2
 */
template <class ScalarT, typename IdxT>
Generator4<ScalarT, IdxT>::Generator4(bus_type* bus, ScalarT P0, ScalarT Q0)
  : ModelEvaluatorImpl<ScalarT, IdxT>(6, 1, 2),
    H_(5.0),
    D_(0.04),
    Xq_(0.85),
    Xd_(1.05),
    Xqp_(0.35),
    Xdp_(0.35),
    Rs_(0.01),
    Tq0p_(1.0), // [s]
    Td0p_(8.0), // [s]
    Ef_(1.45),
    Pm_(1.0),
    omega_s_(1.0),
    omega_b_(2.0*60.0*M_PI),
    omega_up_(omega_s_ + 0.0001),
    omega_lo_(omega_s_ - 0.0001),
    c_(10000.0),
    beta_(2),
    P0_(P0),
    Q0_(Q0),
    bus_(bus)
{
}

template <class ScalarT, typename IdxT>
Generator4<ScalarT, IdxT>::~Generator4()
{
}

/*!
 * @brief This function will be used to allocate sparse Jacobian matrices.
 *
 */
template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::allocate()
{
    //std::cout << "Allocate Generator4..." << std::endl;
    tag_.resize(size_);

    return 0;
}

/**
 * Initialization of the generator model
 *
 */
template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::initialize()
{
    // std::cout << "Initialize Generator4..." << std::endl;

    // Compute initial guess for the generator voltage phase
    const ScalarT delta = atan((Xq_*P0_ - Rs_*Q0_) / (V()*V() + Rs_*P0_ + Xq_*Q0_)) + theta();

    // Compute initial guess for the generator current phase
    const ScalarT phi   = theta() - delta - atan(Q0_/P0_);

    // Compute initial gueses for generator currents and potentials in d-q frame
    const ScalarT Id = std::sqrt(P0_*P0_ + Q0_*Q0_)/V() * sin(phi);
    const ScalarT Iq = std::sqrt(P0_*P0_ + Q0_*Q0_)/V() * cos(phi);
    const ScalarT Ed = V()*sin(theta() - delta) + Rs_*Id + Xqp_*Iq;
    const ScalarT Eq = V()*cos(theta() - delta) + Rs_*Iq - Xdp_*Id;

    y_[0] =  delta;
    y_[1] =  omega_s_;
    y_[2] =  Ed;
    y_[3] =  Eq;
    y_[4] =  Id;
    y_[5] =  Iq;
    yp_[0] = 0.0;
    yp_[1] = 0.0;
    yp_[2] = 0.0;
    yp_[3] = 0.0;
    yp_[4] = 0.0;
    yp_[5] = 0.0;

    // Set control parameter values here.
    Ef_ = Eq - (Xd_ - Xdp_)*Id;                // <~ set to steady state value
    Pm_ = Ed*Id + Eq*Iq + (Xdp_ - Xqp_)*Id*Iq; // <~ set to steady state value

    // Initialize optimization parameters
    param_[0] = Pm_;
    param_up_[0] = 1.5;
    param_lo_[0] = 0.0;

    param_[1] = Ef_;
    param_up_[1] = 1.7;
    param_lo_[1] = 0.0;

    return 0;
}

/**
 * \brief Identify differential variables.
 */
template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::tagDifferentiable()
{
    tag_[0] = true;
    tag_[1] = true;
    tag_[2] = true;
    tag_[3] = true;

    for (IdxT i=4; i < size_; ++i)
    {
        tag_[i] = false;
    }

    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateResidual()
{
    // std::cout << "Evaluate residual for Generator4..." << std::endl;

    f_[0] = -yp_[0] + y_[1]-omega_s_;
    f_[1] = -yp_[1] + omega_s_/(2.0*H_)*( Pm() - y_[3]*y_[5] - y_[2]*y_[4] - (Xdp_ - Xqp_)*y_[4]*y_[5] - D_*(y_[1]-omega_s_) );
    f_[2] = -yp_[2] + 1.0/Tq0p_*( -y_[2] - (Xq_ - Xqp_)*y_[5] );
    f_[3] = -yp_[3] + 1.0/Td0p_*( -y_[3] + (Xd_ - Xdp_)*y_[4] + Ef() );
    f_[4] =  Rs_*y_[4] + Xqp_*y_[5] + V()*sin(theta() - y_[0]) - y_[2];
    f_[5] = -Xdp_*y_[4] + Rs_*y_[5] + V()*cos(theta() - y_[0]) - y_[3];

    // Compute active and reactive load provided by the infinite bus.
    P() = Pg();
    Q() = Qg();
    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateJacobian()
{
    std::cerr << "Evaluate Jacobian for Generator4..." << std::endl;
    std::cerr << "Jacobian evaluation not implemented!" << std::endl;
    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateIntegrand()
{
    // std::cout << "Evaluate Integrand for Generator4..." << std::endl;
    g_[0] = frequencyPenalty(y_[1]);
    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::initializeAdjoint()
{
    //std::cout << "Initialize adjoint for Generator4..." << std::endl;
    for (IdxT i=0; i<size_; ++i)
    {
        yB_[i] = 0.0;
        ypB_[i] = 0.0;
    }
    ypB_[1] = frequencyPenaltyDer(y_[1]);

    return 0;
}

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateAdjointResidual()
{
    // std::cout << "Evaluate adjoint residual for Generator4..." << std::endl;
    fB_[0] = -ypB_[0] + yB_[4]*V()*cos(theta() - y_[0]) - yB_[5]*V()*sin(theta() - y_[0]);
    fB_[1] = -ypB_[1] - yB_[0] + yB_[1]*(omega_s_/(2.0*H_))*D_ + frequencyPenaltyDer(y_[1]);
    fB_[2] = -ypB_[2] + yB_[1]*(omega_s_/(2.0*H_))*y_[4] + yB_[2]/Tq0p_ + yB_[4];
    fB_[3] = -ypB_[3] + yB_[1]*(omega_s_/(2.0*H_))*y_[5] + yB_[3]/Td0p_ + yB_[5];
    fB_[4] = yB_[1]*(omega_s_/(2.0*H_))*(y_[2] + (Xdp_ - Xqp_)*y_[5]) - yB_[3]*(Xd_ - Xdp_)/Td0p_ - yB_[4]*Rs_  + yB_[5]*Xdp_;
    fB_[5] = yB_[1]*(omega_s_/(2.0*H_))*(y_[3] + (Xdp_ - Xqp_)*y_[4]) + yB_[2]*(Xq_ - Xqp_)/Tq0p_ - yB_[4]*Xqp_ - yB_[5]*Rs_;
    return 0;
}

// template <class ScalarT, typename IdxT>
// int Generator4<ScalarT, IdxT>::evaluateAdjointJacobian()
// {
//     std::cout << "Evaluate adjoint Jacobian for Generator4..." << std::endl;
//     std::cout << "Adjoint Jacobian evaluation not implemented!" << std::endl;
//     return 0;
// }

template <class ScalarT, typename IdxT>
int Generator4<ScalarT, IdxT>::evaluateAdjointIntegrand()
{
    // std::cout << "Evaluate adjoint Integrand for Generator4..." << std::endl;
    gB_[0] = -omega_s_/(2.0*H_) * yB_[1];
    gB_[1] = -1.0/Td0p_ * yB_[3];

    return 0;
}


//
// Private functions
//

/**
 * Generator active power Pg.
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::Pg()
{
    return y_[5]*V()*cos(theta() - y_[0]) + y_[4]*V()*sin(theta() - y_[0]);
}

/**
 * Generator reactive power Qg.
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::Qg()
{
    return y_[5]*V()*sin(theta() - y_[0]) - y_[4]*V()*cos(theta() - y_[0]);
}

/**
 * Frequency penalty is used as the objective function for the generator model.
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::frequencyPenalty(ScalarT omega)
{
    return c_ * pow(std::max(0.0, std::max(omega - omega_up_, omega_lo_ - omega)), beta_);
}

/**
 * Derivative of frequency penalty cannot be written in terms of min/max functions.
 * Need to expand conditional statements instead.
 */
template <class ScalarT, typename IdxT>
ScalarT Generator4<ScalarT, IdxT>::frequencyPenaltyDer(ScalarT omega)
{
    if (omega > omega_up_)
    {
        return beta_ * c_ * pow(omega - omega_up_, beta_ - 1.0);
    }
    else if (omega < omega_lo_)
    {
        return beta_ * c_ * pow(omega - omega_lo_, beta_ - 1.0);
    }
    else
    {
        return 0.0;
    }
}



// Available template instantiations
template class Generator4<double, long int>;
template class Generator4<double, size_t>;


} // namespace ModelLib

