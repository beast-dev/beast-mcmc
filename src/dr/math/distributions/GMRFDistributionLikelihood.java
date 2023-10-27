/*
 * GMRFDistributionLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions;

import dr.evomodel.coalescent.OldGMRFSkyrideLikelihood;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class GMRFDistributionLikelihood extends AbstractModelLikelihood implements MultivariateDistribution {

    private Parameter precisionParameter;

    private Parameter lambdaParameter;

    private Parameter gridParameter;

    private int fieldLength;
    private Parameter popSizeParameter;
    private SymmTridiagMatrix weightMatrix;

    private static final double LOG_TWO_TIMES_PI = OldGMRFSkyrideLikelihood.LOG_TWO_TIMES_PI;


    public GMRFDistributionLikelihood(String name,
                                      Parameter precisionParameter,
                                      Parameter lambdaParameter,
                                      Parameter gridParameter,
                                      Parameter popSizeParameter) {
        super(name);
        this.precisionParameter = precisionParameter;
        this.lambdaParameter = lambdaParameter;
        this.gridParameter = gridParameter;
        this.fieldLength = popSizeParameter.getDimension();
        this.popSizeParameter = popSizeParameter;
        addVariable(precisionParameter);
        addVariable(lambdaParameter);
        addVariable(gridParameter);
    }


    @Override
    public double logPdf(double[] x) {
        DenseVector diagonal1 = new DenseVector(x.length);
        DenseVector currentGamma = new DenseVector(x);

        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
        currentQ.mult(currentGamma, diagonal1);

        double currentLike = 0;
        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
        if (lambdaParameter.getParameterValue(0) == 1) {
            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
        } else {
            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
        }

        return currentLike;
    }

    private void setupGMRFWeights() {

        //Set up the weight Matrix
        double[] offdiag = new double[fieldLength - 1];
        double[] diag = new double[fieldLength];

        //First set up the offdiagonal entries;

        for (int i = 0; i < fieldLength - 1; i++) {
            offdiag[i] = -1.0;
        }

        //Then set up the diagonal entries;
        for (int i = 1; i < fieldLength - 1; i++)
            diag[i] = -(offdiag[i] + offdiag[i - 1]);

        //Take care of the endpoints
        diag[0] = -offdiag[0];
        diag[fieldLength - 1] = -offdiag[fieldLength - 2];

        weightMatrix = new SymmTridiagMatrix(diag, offdiag);
    }

    private SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda) {
        if (lambda == 1)
            return getScaledWeightMatrix(precision);

        SymmTridiagMatrix a = weightMatrix.copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
            a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
        }

        a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
        return a;
    }

    private SymmTridiagMatrix getScaledWeightMatrix(double precision) {
        SymmTridiagMatrix a = weightMatrix.copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
        return a;
    }


    private DenseVector getMeanAdjustedGamma(double[] x) {
        DenseVector currentGamma = new DenseVector(x);
        return currentGamma;
    }

    @Override
    public double[][] getScaleMatrix() {
        return new double[0][];
    }

    @Override
    public double[] getMean() {
        return new double[0];
    }

    @Override
    public String getType() {
        return null;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
    public double getLogLikelihood() {
        return 0;
    }

    @Override
    public void makeDirty() {

    }
}
