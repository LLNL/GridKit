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
#include <iomanip>

#include "ComponentLib/Generator2/Generator2.hpp"
#include "Solver/Dynamic/Ida.hpp"

#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>



int main()
{
    using namespace ModelLib;
    using namespace AnalysisManager::Sundials;

    ModelEvaluator<double, double, size_t>* model = new Generator2<double, double, size_t>();
    Ida<double, double, size_t>* idas = new Ida<double, double, size_t>(model);
    
    model->allocate();

    double t_init  = 0.0;  
    double t_final = 50.0;

    // setup simulation
    idas->configureSimulation();
    idas->configureAdjoint();
    idas->getDefaultInitialCondition();
    idas->initializeSimulation(t_init);
    idas->configureQuadrature();
    idas->initializeQuadrature();
    
    double t_fault = 0.1;
    double t_clear = 0.2;
    double g1, g2, eps = 1e-5;
    idas->runSimulation(t_fault);
    // create initial condition after a fault
    {
        Generator2<double, double, size_t>* gen2 = dynamic_cast<Generator2<double, double, size_t>*>(model);
        gen2->shortCircuit();
        idas->runSimulation(t_clear, 2);
        gen2->restore();
        idas->saveInitialCondition();
    }
    
    const double* Q = idas->getIntegral();
    // run simulation
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    idas->initializeQuadrature();
    idas->runSimulationQuadrature(t_final, 10);
    idas->printFinalStats();
    g1 = Q[0];
    std::cout << "Q = " << std::setprecision(15) << std::scientific << g1 << "\n";
    
    // run simulation with perturbed parameter
    model->param()[0] += eps;
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    idas->initializeQuadrature();
    idas->runSimulationQuadrature(t_final);
    idas->printFinalStats();
    g2 = Q[0];
    std::cout << "Q = " << g1 << "\n";
    
    // restore parameter to original value
    model->param()[0] -= eps;
    
    // initialize adjoint and run forward simulation
    idas->initializeAdjoint();
    idas->getSavedInitialCondition();
    idas->initializeSimulation(t_init);
    idas->initializeQuadrature();
    idas->runForwardSimulation(t_final);
    idas->printFinalStats();
    
    // Evaluate dG/dp numerically
    std::cout << "dG/dp (numerical) = " << (g2 - g1)/eps << "\n";
    
    idas->initializeBackwardSimulation(t_final);

    // Backward simulation evaluates dG/dp from adjoint
    idas->runBackwardSimulation(t_init);
    double* neg_dGdp = idas->getAdjointIntegral();
    std::cout << "dG/dp (adjoint)   = " << -neg_dGdp[0]      << "\n";

    idas->deleteSimulation();
    delete idas;
    delete model;

    return 0;
}
