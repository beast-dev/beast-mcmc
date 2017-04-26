/*
 * GMRFSkyrideFixedEffectsGibbsOperator.java
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

package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodelxml.coalescent.operators.GMRFSkyrideFixedEffectsGibbsOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.MultivariateDistribution;
import no.uib.cipr.matrix.*;

/**
 * A Gibbs operator to update the population size parameters under a Gaussian Markov random field prior
 *
 * @author Erik Bloomquist
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineFixedEffectsGibbsOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */

public class GMRFSkyrideFixedEffectsGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private GMRFSkyrideLikelihood gmrfLikelihood;

    private DenseVector mean;
    private DenseMatrix precision;

    private int fieldLength;
    private int dim;

    public GMRFSkyrideFixedEffectsGibbsOperator(Parameter param,
                                                GMRFSkyrideLikelihood gmrfLikelihood, MultivariateDistribution prior, double weight) {
        this.gmrfLikelihood = gmrfLikelihood;
        mean = new DenseVector(prior.getMean());
        precision = new DenseMatrix(prior.getScaleMatrix());

        this.fieldLength = gmrfLikelihood.getPopSizeParameter().getDimension();
        this.dim = param.getDimension();

        this.gmrfLikelihood = gmrfLikelihood;

        setWeight(weight);
    }

    public double doOperation() {

        DenseMatrix X = new DenseMatrix(gmrfLikelihood.getDesignMatrix().getParameterAsMatrix());
        SymmTridiagMatrix Q = gmrfLikelihood.getScaledWeightMatrix(gmrfLikelihood.getPrecisionParameter().getParameterValue(0),
                gmrfLikelihood.getLambdaParameter().getParameterValue(0));
        DenseVector gamma = new DenseVector(gmrfLikelihood.getPopSizeParameter().getParameterValues());

        Parameter.Abstract beta = (Parameter.Abstract) gmrfLikelihood.getBetaParameter();

        //Set up the Vectors and matricies for the gibbs step
        DenseMatrix gibbsPrecision = precision.copy();
        UpperSPDDenseMatrix gibbsVariance;
        DenseVector gibbsMean = new DenseVector(dim);
        DenseMatrix workingMatrix = new DenseMatrix(dim, fieldLength);
        DenseVector workingVector = new DenseVector(dim);

        //Get the correct forms
        X.transAmultAdd(Q, workingMatrix);
        workingMatrix.multAdd(X, gibbsPrecision);

        precision.mult(mean, workingVector);
        workingMatrix.multAdd(gamma, workingVector);

        workingMatrix = Matrices.identity(dim);

        gibbsPrecision.solve(Matrices.identity(dim), workingMatrix);
        gibbsVariance = new UpperSPDDenseMatrix(workingMatrix);
        gibbsVariance.mult(workingVector, gibbsMean);

        //Propose a new value for beta
        DenseVector betaNew = GMRFSkyrideBlockUpdateOperator.getMultiNormal(gibbsMean, gibbsVariance);

        for (int i = 0; i < dim; i++) {
            beta.setParameterValueQuietly(i, betaNew.get(i));
        }

        beta.fireParameterChangedEvent();

        return 0;
    }

    public int getStepCount() {
        return 0;
    }


    public String getPerformanceSuggestion() {
        return null;
    }

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return GMRFSkyrideFixedEffectsGibbsOperatorParser.GMRF_GIBBS_OPERATOR;
    }

}
