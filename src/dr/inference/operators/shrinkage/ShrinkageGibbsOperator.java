/*
 * ShrinkageGibbsOperator.java
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

import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.InverseGammaDistribution;

public class ShrinkageGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{
    IndependentInverseGammaDistributionModel localPrior;
    IndependentInverseGammaDistributionModel globalPrior;
    Parameter data;
    double pathParameter;
    final boolean usePathParameter = true;

    public ShrinkageGibbsOperator(double weight, IndependentInverseGammaDistributionModel localPrior, IndependentInverseGammaDistributionModel globalPrior, Parameter data){
        setWeight(weight);
        this.localPrior = localPrior;
        this.globalPrior = globalPrior;
        this.data = data;
        this.pathParameter = 1;
    }


    public int getStepCount() {
        return 0;
    }

    @Override
    public String getOperatorName() {
        return "ShrinkageGibbsOperator";
    }

    @Override
    public double doOperation() {
        Parameter local = localPrior.getData();
        Parameter localShape = localPrior.getShape();
        Parameter localAugmented = localPrior.getScale();

        Parameter global = globalPrior.getData();
        Parameter globalShape = globalPrior.getShape();
        Parameter globalAugmented = globalPrior.getScale();

        for (int i = 0; i < local.getDimension(); i++) {
            double scale = localAugmented.getParameterValue(i) + data.getParameterValue(i) * data.getParameterValue(i) / (2 * global.getParameterValue(0));
            double draw = InverseGammaDistribution.nextInverseGamma(localShape.getParameterValue(i) + .5, scale);
            local.setParameterValueQuietly(i, draw);
        }
        local.fireParameterChangedEvent();

        double shape = local.getDimension() / 2 + globalShape.getParameterValue(0);
        double scale = globalAugmented.getParameterValue(0);
        for (int i = 0; i < local.getDimension(); i++) {
            scale += .5 * (data.getParameterValue(i) * data.getParameterValue(i)) / local.getParameterValue(i);
        }
        double draw = InverseGammaDistribution.nextInverseGamma(shape, scale);
        global.setParameterValue(0, draw);


        return 0;
    }


}
