/*
 * DiagonalAttenuationGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treedatalikelihood.hmc;

import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public class DiagonalAttenuationGradient extends AbstractDiffusionGradient implements Reportable {

    private final int dim;
    private final BranchSpecificGradient branchSpecificGradient;

    private final DiagonalMatrix attenuation;

    public DiagonalAttenuationGradient(BranchSpecificGradient branchSpecificGradient,
                                       Likelihood likelihood,
                                       MatrixParameterInterface parameter) {
        super(likelihood, 0, Double.POSITIVE_INFINITY);

        assert (parameter instanceof DiagonalMatrix)
                : "DiagonalAttenuationGradient can only be applied to a DiagonalMatrix.";

        this.attenuation = (DiagonalMatrix) parameter;

        this.branchSpecificGradient = branchSpecificGradient;
        this.dim = parameter.getColumnDimension();

    }

    @Override
    public Parameter getParameter() {
        return attenuation.getDiagonalParameter();
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public Parameter getRawParameter() {
        return attenuation;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] gradient = branchSpecificGradient.getGradientLogDensity();
        return getGradientLogDensity(gradient);
    }

    public double[] getGradientLogDensity(double[] gradient) {
        return extractDiagonalGradient(gradient);
    }

    private double[] extractDiagonalGradient(double[] gradient) {
        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = gradient[offset + i];
        }
        return result;
    }

    @Override
    public String getReport() {
        return "attenuationGradient." + attenuation.getParameterName() + "\n" +
                super.getReport();
    }
}