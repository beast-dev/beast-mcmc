/*
 * ModelSpecificPseudoPriorLikelihood.java
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
 */
public class ModelSpecificPseudoPriorLikelihood extends DistributionLikelihood{

    private int[] models;
    protected Distribution prior;
    protected Distribution pseudoPrior;
    private Parameter modelIndicator;

    public ModelSpecificPseudoPriorLikelihood(
            Distribution prior,
            Distribution pseudoPrior,
            Parameter modelIndicator,
            int[] models){
        super(prior);
        this.prior = prior;
        this.pseudoPrior = pseudoPrior;
        this.models = models;
        this.modelIndicator = modelIndicator;
    }


   /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {
        boolean inModel = false;
        int modelCode = (int)modelIndicator.getParameterValue(0);
        for(int i = 0; i < models.length; i++){
            if(models[i] == modelCode){
                inModel = true;
                break;
            }
        }

        if(inModel){
            distribution = prior;
        }else{
            distribution = pseudoPrior;
        }
        double logL = super.calculateLogLikelihood();
        return logL;
    }

}
