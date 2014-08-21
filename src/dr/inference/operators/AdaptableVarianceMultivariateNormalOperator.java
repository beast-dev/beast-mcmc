/*
 * AdaptableVarianceMultivariateNormalOperator.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators;

import java.util.ArrayList;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.util.Transform;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Guy Baele
 * @author Marc A. Suchard
 */
public class AdaptableVarianceMultivariateNormalOperator extends AbstractCoercableOperator {

    public static final String AVMVN_OPERATOR = "adaptableVarianceMultivariateNormalOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String BETA = "beta";
    public static final String INITIAL = "initial";
    public static final String BURNIN = "burnin";
    public static final String UPDATE_EVERY = "updateEvery";
    public static final String FORM_XTX = "formXtXInverse";
    public static final String COEFFICIENT = "coefficient";

    public static final boolean DEBUG = false;

    private double scaleFactor;
    private double beta;
    private int iterations, updates, initial, burnin, every;
    private final Parameter parameter;
    private final Transform[] transformations;
    private final int dim;
    // private final double constantFactor;
    private double[] oldMeans, newMeans;

    final double[][] matrix;
    private double[][] empirical;
    private double[][] cholesky;

    // temporary storage, allocated once.
    private double[] epsilon;
    private double[][] proposal;

    public AdaptableVarianceMultivariateNormalOperator(Parameter parameter, Transform[] transformations, double scaleFactor, double[][] inMatrix,
            double weight, double beta, int initial, int burnin, int every, CoercionMode mode, boolean isVarianceMatrix) {

        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        this.transformations = transformations;
        this.beta = beta;
        this.iterations = 0;
        this.updates = 0;
        setWeight(weight);
        dim = parameter.getDimension();
        // constantFactor = Math.pow(2.38, 2) / ((double) dim); // not necessary because scaleFactor is auto-tuned
        this.initial = initial;
        this.burnin = burnin;
        this.every = every;
        this.empirical = new double[dim][dim];
        this.oldMeans = new double[dim];
        this.newMeans = new double[dim];

        this.epsilon = new double[dim];
        this.proposal = new double[dim][dim];

        SingularValueDecomposition svd = new SingularValueDecomposition(new DenseDoubleMatrix2D(inMatrix));
        if (inMatrix[0].length != svd.rank()) {
            throw new RuntimeException("Variance matrix in AdaptableVarianceMultivariateNormalOperator is not of full rank");
        }

        if (isVarianceMatrix) {
            matrix = inMatrix;
        } else {
            matrix = formXtXInverse(inMatrix);
        }

        /*System.err.println("matrix initialization: ");
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                System.err.print(matrix[i][j] + " ");
            }
            System.err.println();
        }*/

        try {
            cholesky = (new CholeskyDecomposition(matrix)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in AdaptableVarianceMultivariateNormalOperator");
        }
    }

    public AdaptableVarianceMultivariateNormalOperator(Parameter parameter, Transform[] transformations, double scaleFactor,
            MatrixParameter varMatrix, double weight, double beta, int initial, int burnin, int every, CoercionMode mode, boolean isVariance) {
        this(parameter, transformations, scaleFactor, varMatrix.getParameterAsMatrix(), weight, beta, initial, burnin, every, mode, isVariance);
    }

    private double[][] formXtXInverse(double[][] X) {
        int N = X.length;
        int P = X[0].length;

        double[][] matrix = new double[P][P];
        for (int i = 0; i < P; i++) {
            for (int j = 0; j < P; j++) {
                int total = 0;
                for (int k = 0; k < N; k++) {
                    total += X[k][i] * X[k][j];
                }
                matrix[i][j] = total;
            }
        }

        // Take inverse
        matrix = new SymmetricMatrix(matrix).inverse().toComponents();
        return matrix;
    }

    //act as if population mean is known
    private double calculateCovariance(int number, double currentMatrixEntry, double[] values, int firstIndex, int secondIndex) {

        // number will always be > 1 here
        /*double result = currentMatrixEntry * (number - 1);
        result += (values[firstIndex] * values[secondIndex]);
        result += ((number - 1) * oldMeans[firstIndex] * oldMeans[secondIndex] - number * newMeans[firstIndex] * newMeans[secondIndex]);
        result /= ((double) number);*/

        double result = currentMatrixEntry * (number - 2);
        result += (values[firstIndex] * values[secondIndex]);
        result += ((number - 1) * oldMeans[firstIndex] * oldMeans[secondIndex] - number * newMeans[firstIndex] * newMeans[secondIndex]);
        result /= ((double)(number - 1));

        return result;

    }

    public double doOperation() throws OperatorFailedException {

        iterations++;

        if (DEBUG) {
            System.err.println("\nAVMVN Iteration: " + iterations);
            System.err.println("Using AdaptableVarianceMultivariateNormalOperator: " + iterations + " for " + parameter.getParameterName());
            System.err.println("Old parameter values:");
            for (int i = 0; i < dim; i++) {
                System.err.println(parameter.getParameterValue(i));
            }
        }

        double[] x = parameter.getParameterValues();

        //transform to the appropriate scale
        double[] transformedX = new double[dim];
        for (int i = 0; i < dim; i++) {
            transformedX[i] = transformations[i].transform(x[i]);
        }

        if (DEBUG) {
            System.err.println("Old transformed parameter values:");
            for (int i = 0; i < dim; i++) {
                System.err.println(transformedX[i]);
            }
        }

        //store MH-ratio in logq
        double logJacobian = 0.0;

        //change this: make a rule for when iterations == burnin
        if (iterations > 1 && iterations > burnin) {

            if (DEBUG) {
                System.err.println("  AVMVN iterations > burnin");
            }

            if (iterations > (burnin+1)) {

                if (iterations % every == 0) {

                    updates++;

                    if (DEBUG) {
                        System.err.println("updates = " + updates);
                    }

                    //first recalculate the means using recursion
                    for (int i = 0; i < dim; i++) {
                        newMeans[i] = ((oldMeans[i] * (updates - 1)) + transformedX[i]) / updates;
                    }

                    if (updates > 1) {
                        //here we can simply use the double[][] matrix
                        for (int i = 0; i < dim; i++) {
                            for (int j = i; j < dim; j++) {
                                empirical[i][j] = calculateCovariance(updates, empirical[i][j], transformedX, i, j);
                                empirical[j][i] = empirical[i][j];
                            }
                        }
                    }

                    if (DEBUG) {
                        System.err.println("Old means:");
                        for (int i = 0; i < dim; i++) {
                            System.err.println(oldMeans[i]);
                        }
                        System.err.println("New means:");
                        for (int i = 0; i < dim; i++) {
                            System.err.println(newMeans[i]);
                        }
                        System.err.println("Empirical covariance matrix:");
                        for (int i = 0; i < dim; i++) {
                            for (int j = 0; j < dim; j++) {
                                System.err.print(empirical[i][j] + " ");
                            }
                            System.err.println();
                        }
                    }

                }

                /*if (iterations == 17) {
                    System.exit(0);
                }*/

            } else if (iterations == (burnin+1)) {

                //updates++;

                //i.e. iterations == burnin+1, i.e. first sample for C_t
                //this will not be reached when burnin is set to 0
                for (int i = 0; i < dim; i++) {
                    //oldMeans[i] = transformedX[i];
                    //newMeans[i] = transformedX[i];
                    oldMeans[i] = 0.0;
                    newMeans[i] = 0.0;
                }

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        empirical[i][j] = 0.0;
                    }
                }

            }

        } else if (iterations == 1) {

            if (DEBUG) {
                System.err.println("  iterations == 1");
            }
            //System.err.println("Iteration: " + iterations);

            //iterations == 1
            for (int i = 0; i < dim; i++) {
                //oldMeans[i] = transformedX[i];
                //newMeans[i] = transformedX[i];
                oldMeans[i] = 0.0;
                newMeans[i] = 0.0;
            }

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    empirical[i][j] = 0.0;
                    proposal[i][j] = matrix[i][j];
                }
            }

        }

        for (int i = 0; i < dim; i++) {
            epsilon[i] = scaleFactor * MathUtils.nextGaussian();
        }

        if (iterations > initial) {

            if (DEBUG) {
                System.err.println("  iterations > initial");
            }

            if (iterations % every == 0) {
                // TODO: For speed, it may not be necessary to update decomposition each and every iteration
                //double start = System.nanoTime();

                // double[][] proposal = new double[dim][dim];
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) { // symmetric matrix
                        proposal[j][i] = proposal[i][j] = (1 - beta) * // constantFactor *  /* auto-tuning using scaleFactor */
                                empirical[i][j] + beta * matrix[i][j];
                    }
                }

                // not necessary for first test phase, but will need to be performed when covariance matrix is being updated
                try {
                    cholesky = (new CholeskyDecomposition(proposal)).getL();
                } catch (IllegalDimension illegalDimension) {
                    throw new RuntimeException("Unable to decompose matrix in AdaptableVarianceMultivariateNormalOperator");
                }

                //double end = System.nanoTime();
                //double baseResult = end - start;
                //System.err.println("Cholesky decomposition took: " + baseResult);
            }

        }

        if (DEBUG) {
            System.err.println("  Drawing new values");
        }

        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                transformedX[i] += cholesky[j][i] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
            if (MULTI) {
                parameter.setParameterValueQuietly(i, transformations[i].inverse(transformedX[i]));
            } else {
                parameter.setParameterValue(i, transformations[i].inverse(transformedX[i]));
            }

            //this should be correct
            //logJacobian += transformations[i].getLogJacobian(parameter.getParameterValue(i)) - transformations[i].getLogJacobian(x[i]);
            logJacobian += transformations[i].getLogJacobian(x[i]) - transformations[i].getLogJacobian(parameter.getParameterValue(i));

        }

        if (DEBUG) {
            for (int i = 0; i < dim; i++) {
                System.err.println(x[i] + " -> " + parameter.getValue(i));
            }
        }

        if (MULTI) {
            parameter.fireParameterChangedEvent(); // Signal once.
        }

        if (iterations % every == 0) {
            if (DEBUG) {
                System.err.println("  Copying means");
            }
            //copy new means to old means for next update iteration
            //System.arraycopy(newMeans, 0, oldMeans, 0, dim);
            double[] tmp = oldMeans;
            oldMeans = newMeans;
            newMeans = tmp; // faster to swap pointers
        }

        //System.err.println("scale factor: " + scaleFactor);
        /*System.err.println("New parameter values:");
        for (int i = 0; i < dim; i++) {
                System.err.println(parameter.getParameterValue(i));
        }*/
        //System.err.println("log(Jacobian): " + logJacobian);

        //return 0.0;
        return logJacobian;

    }

    public String toString() {
        return AVMVN_OPERATOR + "(" + parameter.getParameterName() + ")";
    }

    public static final boolean MULTI = true;

    //Methods needed when using TwoPhaseOperator(Parser)
    public Parameter getParameter() {
        return this.parameter;
    }

    public void provideSamples(ArrayList<ArrayList<Double>> parameterSamples) {
        if (DEBUG) {
            System.err.println("AVMVN operator parameter length: " + parameter.getDimension());
            System.err.println("provideSamples argument length: " + parameterSamples.size());
        }
        if (parameter.getDimension() != parameterSamples.size()) {
            throw new RuntimeException("Dimension mismatch in AVMVN Operator: inconsistent parameter dimensions");
        } else {
            int lowestNumberOfSamples = parameterSamples.get(0).size();
            for (int i = 0; i < parameterSamples.size(); i++) {
                if (parameterSamples.get(i).size() < lowestNumberOfSamples) {
                    lowestNumberOfSamples = parameterSamples.get(i).size();
                }
            }
            if (DEBUG) {
                System.err.println("lowest number of samples: " + lowestNumberOfSamples);
            }
            //set number of iterations of AVMVN operator
            this.iterations = lowestNumberOfSamples;
            this.updates = lowestNumberOfSamples;
            this.beta = 0.0;
            //set means based on provided samples, but take into account transformation(s)
            for (int i = 0; i < parameterSamples.size(); i++) {
                for (int j = 0; j < lowestNumberOfSamples; j++) {
                    newMeans[i] += transformations[i].transform(parameterSamples.get(i).get(j));
                    //parameterSamples.get(i).get(j);
                }
                newMeans[i] /= (double)lowestNumberOfSamples;
            }
            if (DEBUG) {
                System.err.println();
                for (int i = 0; i < parameterSamples.size(); i++) {
                    System.err.println("Mean " + i + ": " + newMeans[i]);
                }
            }
            //set covariance matrix based on provided samples, but take into account transformation(s)
            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    for (int k = 0; k < lowestNumberOfSamples; k++) {
                        empirical[i][j] += transformations[i].transform(parameterSamples.get(i).get(k))*transformations[i].transform(parameterSamples.get(j).get(k));
                    }
                    empirical[i][j] /= (double)lowestNumberOfSamples;
                    empirical[i][j] -= newMeans[i]*newMeans[j];
                    empirical[j][i] = empirical[i][j];
                }
            }
            if (DEBUG) {
                System.err.println();
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        System.err.print(empirical[i][j] + "  ");
                    }
                    System.err.println();
                }
            }
        }
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "adaptableVarianceMultivariateNormal(" + parameter.getParameterName() + ")";
    }

    public double getCoercableParameter() {
        return Math.log(scaleFactor);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = Math.exp(value);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return AVMVN_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            if (DEBUG) {
                System.err.println("Parsing AdaptableVarianceMultivariateNormalOperator.");
            }

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double beta = xo.getDoubleAttribute(BETA);
            int initial = xo.getIntegerAttribute(INITIAL);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
            double coefficient = xo.getDoubleAttribute(COEFFICIENT);

            int burnin = 0;
            int every = 1;
            if (xo.hasAttribute(BURNIN)) {
                burnin = xo.getIntegerAttribute(BURNIN);
            }
            if (burnin > initial || burnin < 0) {
                throw new XMLParseException("burnin must be smaller than the initial period");
            }

            if (xo.hasAttribute(UPDATE_EVERY)) {
                every = xo.getIntegerAttribute(UPDATE_EVERY);
            }

            if (every <= 0) {
                throw new XMLParseException("covariance matrix needs to be updated at least every single iteration");
            }

            if (scaleFactor <= 0.0) {
                throw new XMLParseException("scaleFactor must be greater than 0.0");
            }

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            boolean formXtXInverse = xo.getAttribute(FORM_XTX, false);

            //varMatrix needs to be initialized
            int dim = parameter.getDimension();
            System.err.println("Dimension: " + dim);

            if (initial <= 2 * dim) {
                initial = 2 * dim;
            }

            Parameter[] init = new Parameter[dim];
            for (int i = 0; i < dim; i++) {
                init[i] = new Parameter.Default(dim, 0.0);
            }
            for (int i = 0; i < dim; i++) {
                init[i].setParameterValue(i, Math.pow(coefficient, 2) / ((double) dim));
            }
            MatrixParameter varMatrix = new MatrixParameter(null, init);

            /*double[][] test = varMatrix.getParameterAsMatrix();
            for (int i = 0; i < dim; i++) {
            	for (int j = 0; j < dim; j++) {
            		System.err.print(test[i][j] + " ");
            	}
            	System.err.println();
            }*/

            Transform[] transformations = new Transform[dim];
            for (int i = 0; i < dim; i++) {
                transformations[i] = Transform.NONE;
            }

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Transform.ParsedTransform) {
                    Transform.ParsedTransform thisObject = (Transform.ParsedTransform) child;

                    if (DEBUG) {
                        System.err.println("Transformations:");
                    }
                    for (int j = thisObject.start; j < thisObject.end; ++j) {
                        transformations[j] = thisObject.transform;
                        if (DEBUG) {
                            System.err.print(transformations[j].getTransformName() + " ");
                        }
                    }
                    if (DEBUG) {
                        System.err.println();
                    }
                }
            }

            // Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

            if (!formXtXInverse) {
                if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
                    throw new XMLParseException("The variance matrix is not square");
            }

            if (varMatrix.getColumnDimension() != parameter.getDimension())
                throw new XMLParseException("The parameter and variance matrix have differing dimensions");

            /*java.util.logging.Logger.getLogger("dr.inference").info("\nCreating the adaptable variance multivariate normal operator:" +
					"\n beta = " + beta + "\n initial = " + initial + "\n burnin = " + burnin + "\n every = " + every +
					"\n If you use this operator, please cite: " + 
			"   Guy Baele, Philippe Lemey, Marc A. Suchard. 2014. In preparation.");*/

            return new AdaptableVarianceMultivariateNormalOperator(parameter, transformations, scaleFactor, varMatrix, weight, beta, initial, burnin, every, mode, !formXtXInverse);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an adaptable variance multivariate normal operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newDoubleRule(BETA),
                AttributeRule.newDoubleRule(COEFFICIENT),
                AttributeRule.newIntegerRule(INITIAL),
                AttributeRule.newIntegerRule(BURNIN, true),
                AttributeRule.newIntegerRule(UPDATE_EVERY, true),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(FORM_XTX, true),
                new ElementRule(Parameter.class),
                new ElementRule(Transform.ParsedTransform.class, 0, Integer.MAX_VALUE)
        };

    };
}
