/*
 * CompoundParameterGradient.java
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

package dr.inference.hmc;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * Gradient wrt a CompoundParameter whose blocks (sub-parameters) have their own gradients, concatenated in block
 * order. Unlike {@link CompoundGradient}, the parameter is the given CompoundParameter itself, and it may carry a
 * simplex sum: e.g. allNus, whose blocks are the relative rates of different partitions and which sums to one.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class CompoundParameterGradient implements GradientWrtParameterProvider, Reportable {

    private final CompoundParameter parameter;
    private final List<GradientWrtParameterProvider> blocks;
    private final Likelihood likelihood;
    private final double simplexSum;
    private final Double tolerance;

    /**
     * @param blocks     one gradient per sub-parameter of parameter, in order, each wrt that sub-parameter
     * @param simplexSum the fixed sum of the parameter, or NaN if it is unconstrained
     */
    public CompoundParameterGradient(CompoundParameter parameter, List<GradientWrtParameterProvider> blocks,
                                     double simplexSum, Double tolerance) {
        if (blocks.size() != parameter.getParameterCount()) {
            throw new IllegalArgumentException(parameter.getId() + " has " + parameter.getParameterCount()
                    + " blocks but " + blocks.size() + " gradients were given");
        }
        List<Likelihood> likelihoods = new ArrayList<>();
        for (int i = 0; i < blocks.size(); ++i) {
            GradientWrtParameterProvider block = blocks.get(i);
            if (block.getParameter() != parameter.getParameter(i)
                    || block.getDimension() != parameter.getParameter(i).getDimension()) {
                throw new IllegalArgumentException("Gradient " + i + " is not wrt block " + i + " of "
                        + parameter.getId());
            }
            for (Likelihood member : block.getLikelihood().getLikelihoodSet()) {
                if (!likelihoods.contains(member)) {
                    likelihoods.add(member);
                }
            }
        }
        this.parameter = parameter;
        this.blocks = blocks;
        this.likelihood = likelihoods.size() == 1 ? likelihoods.get(0) : new CompoundLikelihood(likelihoods);
        this.simplexSum = simplexSum;
        this.tolerance = tolerance;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public GradientType getGradientType() {
        return GradientType.combine(blocks);
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        final double[] gradient = new double[getDimension()];
        int offset = 0;
        for (GradientWrtParameterProvider block : blocks) {
            final double[] blockGradient = block.getGradientLogDensity();
            System.arraycopy(blockGradient, 0, gradient, offset, blockGradient.length);
            offset += blockGradient.length;
        }
        return gradient;
    }

    @Override
    public double getSimplexSum() {
        return simplexSum;
    }

    // a block whose likelihood is only smooth on coarse steps sets the step for the whole check
    @Override
    public double getNumericGradientStepSize() {
        double step = GradientWrtParameterProvider.super.getNumericGradientStepSize();
        for (GradientWrtParameterProvider block : blocks) {
            step = Math.max(step, block.getNumericGradientStepSize());
        }
        return step;
    }

    @Override
    public String getReport() {
        return "compoundParameterGradient." + parameter.getId() + "\n"
                + GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, tolerance);
    }
}
