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

#include <ModelEvaluatorImpl.hpp>

namespace ModelLib
{
    /*!
     * @brief Implementation of a second order generator model.
     * 
     */
    template  <class ScalarT, typename T, typename I>
    class Generator2 : public ModelEvaluatorImpl<ScalarT, T, I>
    {
        using ModelEvaluatorImpl<ScalarT, T, I>::size_;
        using ModelEvaluatorImpl<ScalarT, T, I>::nnz_;
        using ModelEvaluatorImpl<ScalarT, T, I>::time_;
        using ModelEvaluatorImpl<ScalarT, T, I>::alpha_;
        using ModelEvaluatorImpl<ScalarT, T, I>::rtol_;
        using ModelEvaluatorImpl<ScalarT, T, I>::atol_;
        using ModelEvaluatorImpl<ScalarT, T, I>::y_;
        using ModelEvaluatorImpl<ScalarT, T, I>::yp_;
        using ModelEvaluatorImpl<ScalarT, T, I>::f_;
        using ModelEvaluatorImpl<ScalarT, T, I>::g_;
        using ModelEvaluatorImpl<ScalarT, T, I>::yB_;
        using ModelEvaluatorImpl<ScalarT, T, I>::ypB_;
        using ModelEvaluatorImpl<ScalarT, T, I>::fB_;
        using ModelEvaluatorImpl<ScalarT, T, I>::gB_;
        using ModelEvaluatorImpl<ScalarT, T, I>::param_;
        using ModelEvaluatorImpl<ScalarT, T, I>::param_up_;
        using ModelEvaluatorImpl<ScalarT, T, I>::param_lo_;

    public:
        Generator2();
        virtual ~Generator2();
        
        int allocate();
        int initialize();
        int tagDifferentiable() {return -1;}
        int evaluateResidual();
        int evaluateJacobian();
        int evaluateIntegrand();
        
        int initializeAdjoint();
        int evaluateAdjointResidual();
        //int evaluateAdjointJacobian();
        int evaluateAdjointIntegrand();
        
        // Generator2 specific public methods
        int shortCircuit();
        int restore();
        
    private:
        inline ScalarT frequencyPenalty(ScalarT omega);
        inline ScalarT frequencyPenaltyDer(ScalarT omega);
        
    private:
        T H_;
        T D_;
        T Pm_;
        T Pmax_; // [Wb]
        T omega_s_;
        T omega_b_;
        T omega_up_;
        T omega_lo_;
        T theta_s_;
        T c_;
        T beta_;
        
    };
    
} // namespace ModelLib