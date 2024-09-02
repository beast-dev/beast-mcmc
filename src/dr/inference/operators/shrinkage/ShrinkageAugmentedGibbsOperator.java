/*
 * ShrinkageAugmentedGibbsOperator.java
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

package dr.inference.operators.shrinkage;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.InverseGammaDistribution;

public class ShrinkageAugmentedGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    Parameter localShape;
    Parameter globalShape;

    Parameter localAugmented;
    Parameter globalAugmented;
    Parameter local;
    Parameter global;

    double localAugmentedShape;
    double localAugmentedRate;
    double globalAugmentedShape;
    double globalAugmentedRate;

    double pathParameter;


    public ShrinkageAugmentedGibbsOperator(double weight, DistributionLikelihood localAugmentedPrior, DistributionLikelihood globalAugmentedPrior,
                                           IndependentInverseGammaDistributionModel localPrior,
                                           IndependentInverseGammaDistributionModel globalPrior){
        setWeight(weight);
        GammaDistribution localAugmentedPriorGamma = (GammaDistribution) localAugmentedPrior.getDistribution();
        GammaDistribution globalAugmentedPriorGamma = (GammaDistribution) globalAugmentedPrior.getDistribution();

        this.localShape = localPrior.getShape();
        this.globalShape = globalPrior.getShape();

        this.local = localPrior.getData();
        this.global = globalPrior.getData();
        this.localAugmented = localPrior.getScale();
        this.globalAugmented = globalPrior.getScale();

        localAugmentedShape = localAugmentedPriorGamma.getShape();
        localAugmentedRate = 1 / localAugmentedPriorGamma.getScale();
        globalAugmentedShape = globalAugmentedPriorGamma.getShape();
        globalAugmentedRate = 1 / globalAugmentedPriorGamma.getScale();
        this.pathParameter = 1;
    }

    public int getStepCount() {
        return 0;
    }

    @Override
    public String getOperatorName() {
        return "ShrinkageAugmentedGibbsOperator";
    }

    @Override
    public double doOperation() {
        for (int i = 0; i < local.getDimension(); i++) {
            double shape = localAugmentedShape +  localShape.getParameterValue(i);
            double rate = localAugmentedRate +  1 / local.getParameterValue(i);
//            double scale = 1 / localAugmentedRate + pathParameter * local.getParameterValue(i);
            localAugmented.setParameterValueQuietly(i, GammaDistribution.nextGamma(shape, 1 / rate));
        }
        localAugmented.fireParameterChangedEvent();

        double shape = globalAugmentedShape +  globalShape.getParameterValue(0);
        double rate = globalAugmentedRate +  1 / global.getParameterValue(0);
//        double scale = 1 / globalAugmentedRate + pathParameter *  global.getParameterValue(0);
        globalAugmented.setParameterValue(0, GammaDistribution.nextGamma(shape, 1 / rate));

        return 0;
    }

}
