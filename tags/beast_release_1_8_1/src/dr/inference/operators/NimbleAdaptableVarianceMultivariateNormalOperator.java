/*
 * NimbleAdaptableVarianceMultivariateNormalOperator.java
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
 */
public class NimbleAdaptableVarianceMultivariateNormalOperator extends AbstractCoercableOperator {

    public static final String NIMBLE_AVMVN_OPERATOR = "nimbleAdaptableVarianceMultivariateNormalOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String UPDATE_EVERY = "updateEvery";
    public static final String FORM_XTX = "formXtXInverse";

    public static final boolean DEBUG = false;

    private double scaleFactor;
    private int iterations, every;
    private final Parameter parameter;
    private final Transform[] transformations;
    private final int dim;
    private double acceptanceRate;

    private double[] statSums;
    private double[][] statProds;
    private int timesAdapted, timesRan;

    final double[][] originalProposalMatrix;
    private double[][] empirical;
    private double[][] cholesky;

    // temporary storage, allocated once.
    private double[] epsilon;
    private double[][] proposal;

    public NimbleAdaptableVarianceMultivariateNormalOperator(Parameter parameter, Transform[] transformations, double scaleFactor, double[][] inMatrix,
            double weight, int every, CoercionMode mode, boolean isVarianceMatrix) {

        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        this.transformations = transformations;
        this.iterations = 0;

        setWeight(weight);

        this.dim = parameter.getDimension();
        if (this.dim ==1 ) {
            this.acceptanceRate = 0.44;
        } else if (this.dim == 2) {
            this.acceptanceRate = 0.35;
        } else if (this.dim == 3) {
            this.acceptanceRate = 0.32;
        } else if (this.dim == 4) {
            this.acceptanceRate = 0.25;
        } else {
            this.acceptanceRate = 0.234;
        }

        this.every = every;

        this.statSums = new double[dim];
        this.statProds = new double[dim][dim];
        this.timesAdapted = 0;
        this.timesRan = 0;

        this.epsilon = new double[dim];
        this.proposal = new double[dim][dim];
        this.empirical = new double[dim][dim];

        SingularValueDecomposition svd = new SingularValueDecomposition(new DenseDoubleMatrix2D(inMatrix));
        if (inMatrix[0].length != svd.rank()) {
            throw new RuntimeException("Variance matrix in AdaptableVarianceMultivariateNormalOperator is not of full rank");
        }

        if (isVarianceMatrix) {
            originalProposalMatrix = inMatrix;
        } else {
            originalProposalMatrix = formXtXInverse(inMatrix);
        }

        if (DEBUG) {
            System.err.println("matrix initialization: ");
            for (int i = 0; i < originalProposalMatrix.length; i++) {
                for (int j = 0; j < originalProposalMatrix.length; j++) {
                    System.err.print(originalProposalMatrix[i][j] + " ");
                }
                System.err.println();
            }
        }

        try {
            cholesky = (new CholeskyDecomposition(originalProposalMatrix)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in NimbleAdaptableVarianceMultivariateNormalOperator");
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                empirical[i][j] = 0.0;
            }
        } 

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                proposal[i][j] = originalProposalMatrix[i][j];
            }
        } 

    }

    public NimbleAdaptableVarianceMultivariateNormalOperator(Parameter parameter, Transform[] transformations, double scaleFactor,
            MatrixParameter varMatrix, double weight, int every, CoercionMode mode, boolean isVariance) {
        this(parameter, transformations, scaleFactor, varMatrix.getParameterAsMatrix(), weight, every, mode, isVariance);
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

    public double doOperation() throws OperatorFailedException {

        iterations++;
        timesRan++;

        if (DEBUG) {
            System.err.println("\nNAVMVN Iteration: " + iterations);

            System.err.println("Using NimbleAdaptableVarianceMultivariateNormalOperator: " + iterations + " for " + parameter.getParameterName());
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

        for (int i = 0; i < dim; i++) {
            epsilon[i] = scaleFactor * MathUtils.nextGaussian();
        }

        for (int i = 0; i < dim; i++) {
            statSums[i] += transformedX[i];
            for (int j = i; j < dim; j++) {
                statProds[i][j] += transformedX[i]*transformedX[j];
                statProds[j][i] = statProds[i][j];
            }
        }

        if (DEBUG) {
            System.err.println("Sliding window approach\nstatSums:");
            for (int i = 0; i < dim; i++) {
                System.err.println(statSums[i]);
            }
            System.err.println("statProds:");
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    System.err.print(statProds[i][j] + "  ");
                }
                System.err.println();
            }
            double[][] tempMatrix = new double[dim][dim];
            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    tempMatrix[i][j] = (statProds[i][j] - statSums[i]*statSums[j]/(double)timesRan)/((double)(timesRan-1.0));
                    tempMatrix[j][i] = tempMatrix[i][j];
                }
            }
            System.err.println("empirCov:");
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    System.err.print(tempMatrix[i][j] + "  ");
                }
                System.err.println();
            }
        }

        if (iterations % every == 0) {

            if (DEBUG) {
                System.err.println("Adapting at iteration: " + iterations);
            }

            timesAdapted++;

            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    empirical[i][j] = (statProds[i][j] - statSums[i]*statSums[j]/(double)timesRan)/((double)(timesRan-1.0));
                    empirical[j][i] = empirical[i][j];
                }
            }

            if (DEBUG) {
                System.err.println("gamma1: " + 1.0/Math.pow((timesAdapted + 3.0),0.8));
            }

            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) { // symmetric matrix
                    proposal[i][j] = proposal[i][j] + (1.0/Math.pow((timesAdapted + 3.0),0.8)) * (empirical[i][j] - proposal[i][j]);
                    proposal[j][i] = proposal[i][j];
                }
            }

            // not necessary for first test phase, but will need to be performed when covariance matrix is being updated
            try {
                cholesky = (new CholeskyDecomposition(proposal)).getL();
            } catch (IllegalDimension illegalDimension) {
                throw new RuntimeException("Unable to decompose matrix in AdaptableVarianceMultivariateNormalOperator");
            }

            //reset sliding window
            timesRan = 0;
            for (int i = 0; i < dim; i++) {
                statSums[i] = 0.0;
                for (int j = 0; j < dim; j++) {
                    statProds[i][j] = 0.0;
                }
            }

        }

        if (DEBUG) {
            System.err.println("Proposal matrix:");
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    System.err.print(proposal[i][j] + "  ");
                }
                System.err.println();
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

        return logJacobian;

    }

    public String toString() {
        return NIMBLE_AVMVN_OPERATOR + "(" + parameter.getParameterName() + ")";
    }

    public static final boolean MULTI = true;

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "nimbleAdaptableVarianceMultivariateNormal(" + parameter.getParameterName() + ")";
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
        //make dependent upon dimension
        return this.acceptanceRate;
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
            return NIMBLE_AVMVN_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            if (DEBUG) {
                System.err.println("Parsing NimbleAdaptableVarianceMultivariateNormalOperator.");
            }

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            int every = 200;
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
            Parameter[] init = new Parameter[dim];
            for (int i = 0; i < dim; i++) {
                init[i] = new Parameter.Default(dim, 0.0);
            }
            for (int i = 0; i < dim; i++) {
                init[i].setParameterValue(i, 1.0/((double)dim));
            }
            MatrixParameter varMatrix = new MatrixParameter(null, init);

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

            return new NimbleAdaptableVarianceMultivariateNormalOperator(parameter, transformations, scaleFactor, varMatrix, weight, every, mode, !formXtXInverse);
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
                AttributeRule.newIntegerRule(UPDATE_EVERY, true),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(FORM_XTX, true),
                new ElementRule(Parameter.class),
                new ElementRule(Transform.ParsedTransform.class, 0, Integer.MAX_VALUE)
        };

    };
}
