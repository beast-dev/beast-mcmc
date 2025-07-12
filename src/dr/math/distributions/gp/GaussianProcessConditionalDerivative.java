/*
 * GaussianProcessPrediction.java
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

package dr.math.distributions.gp;

import dr.inference.model.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static dr.math.distributions.gp.AdditiveGaussianProcessDistribution.BasisDimension;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class GaussianProcessConditionalDerivative extends GaussianProcessPrediction {

    private final Parameter realizedValues;
    private final int dim;
    private final Parameter priorMeanDerivative;


    public GaussianProcessConditionalDerivative(AdditiveGaussianProcessDistribution gp,
                                                Parameter realizedValues,
                                                Parameter priorMeanDerivative) {
        super(gp, realizedValues, extractDesignMatrices(gp));

        this.gp = gp;
        this.dim = realizedValues.getDimension();
        this.priorMeanDerivative = priorMeanDerivative;

        this.realizedValues = realizedValues;
        this.mean = new double[dim];

        priorMeanDerivative.addVariableListener(this);
    }

    public GaussianProcessConditionalDerivative(AdditiveGaussianProcessDistribution gp,
                                                Parameter realizedValues) {
        this(gp, realizedValues, setDefaultPriorMean(realizedValues.getDimension()));
    }

    private static Parameter setDefaultPriorMean(int dimension) {
        Parameter priorMeanDerivative = null;;
        for (int i = 0; i < dimension; i++) {
            priorMeanDerivative.setParameterValue(i, 0.0);
        }
        return priorMeanDerivative;
    }

    private static List<DesignMatrix> extractDesignMatrices(AdditiveGaussianProcessDistribution gp) {
        List<DesignMatrix> designMatrices = new ArrayList<>();
        for (BasisDimension basis : gp.getBases()) {
            designMatrices.add(basis.getDesignMatrix1());
        }
        return designMatrices;
    }

    protected List<BasisDimension> makeBases(boolean doSecondDerivative) {
        List<BasisDimension> result = new ArrayList<>();
        List<BasisDimension> originalBases = gp.getBases();

        for (int i = 0; i < originalBases.size(); ++i) {
            BasisDimension originalBasis = originalBases.get(i);
            GaussianProcessKernel originalKernel = originalBasis.getKernel();
            DesignMatrix designMatrix = originalBasis.getDesignMatrix1();

            GaussianProcessKernel.KernelDerivatives kernel =
                    new GaussianProcessKernel.KernelDerivatives(originalKernel, doSecondDerivative);
            BasisDimension newBasis = new BasisDimension(kernel, designMatrix, designMatrix);

            result.add(newBasis);
        }

        return result;
    }

    protected void computeMean() {
        if (!meanKnown) {
            computeCrossRealized();
            double[] meanOriginal =  gp.getMean();
            mean = Arrays.copyOf(priorMeanDerivative.getParameterValues(), priorMeanDerivative.getDimension());
            for(int i = 0; i < dim; i += 1) {
                for(int j = 0; j < dim; j += 1) {
                    mean[i] += crossRealized.get(i,j) * (realizedValues.getParameterValue(j) - meanOriginal[j]);
                }
            }
            meanKnown = true;
        }
    }

//    @Override
//    public void modelChangedEvent(Model model, Object object, int index) { //TODO this should build again the cross and realized bases!!!
//        if (model == gp) {
//            predictionKnown = false;
//            meanKnown = false;
//            varianceKnown = false;
//            crossRealizedKnown = false;
//            crossGramianKnown = false;
//        } else {
//            throw new IllegalArgumentException("Unknown model");
//        }
//    }

//    @Override
//    public void modelRestored(Model model) { predictionKnown = false; }

//    @Override //TODO
//    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
//        if (variable == realizedValues) {
//            predictionKnown = false;
//            meanKnown = false;
//        } else if (variable instanceof DesignMatrix &&
//                predictiveDesigns.contains((DesignMatrix) variable)) {
//            throw new IllegalArgumentException("Not yet implemented");
//        } else {
//            throw new IllegalArgumentException("Unknown variable");
//        }
//    }
}