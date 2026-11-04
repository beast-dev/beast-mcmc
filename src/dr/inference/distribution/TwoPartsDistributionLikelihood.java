/*
 * TwoPartsDistributionLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.inference.distribution;

import dr.math.distributions.Distribution;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * This class facilitates the pseudo-prior for model averaging by bayesian stochastic sampling.
 */
public class TwoPartsDistributionLikelihood extends DistributionLikelihood{
    public static final int PRESENT = 1;
    public static final int ABSENT = 0;
    protected Distribution prior;
    protected Distribution pseudoPrior;

    protected Parameter bitVector;
    protected int paramIndex;

    public TwoPartsDistributionLikelihood(
            Distribution prior,
            Distribution pseudoPrior,
            Parameter bitVector,
            int paramIndex){
        super(prior);
        this.prior = distribution;
        this.pseudoPrior = pseudoPrior;
        this.bitVector = bitVector;
        this.paramIndex = paramIndex;
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        int paramStatus = (int)bitVector.getParameterValue(paramIndex);
        //System.out.println(paramStatus);
        if(paramStatus == PRESENT){
            distribution = prior;
        }else if(paramStatus == ABSENT){
            distribution = pseudoPrior;
        }
        double logL = super.calculateLogLikelihood();
        //System.out.println(paramStatus+" "+logL);
        return logL;
    }
}
