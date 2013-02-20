/*
 * AdaptiveMetropolisOperator.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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
import dr.xml.*;


/**
 * @author Guy Baele
 */
public class AdaptiveMetropolisOperator extends AbstractCoercableOperator {

    public static final String AM_OPERATOR = "adaptiveMetropolisOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String BETA = "beta";
    public static final String FORM_XTX = "formXtXInverse";

    private double scaleFactor;
    private double beta;
    private int iterations, cutoff;
    private final Parameter parameter;
    private final int dim;
    private double[] oldMeans, newMeans;

    final double[][] matrix;
    private double[][] empirical;
    private double[][] cholesky;

    public AdaptiveMetropolisOperator(Parameter parameter, double scaleFactor, double[][] inMatrix, double weight,
                                      double beta, CoercionMode mode, boolean isVarianceMatrix) {

        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        this.beta = beta;
        this.iterations = 0;
        setWeight(weight);
        dim = parameter.getDimension();
        this.cutoff = 2*dim;
        this.empirical = new double[dim][dim];
        this.oldMeans = new double[dim];
        this.newMeans = new double[dim];

        SingularValueDecomposition svd = new SingularValueDecomposition(new DenseDoubleMatrix2D(inMatrix));
        if (inMatrix[0].length != svd.rank()) {
            throw new RuntimeException("Variance matrix in mvnOperator is not of full rank");
        }

        if (isVarianceMatrix) {
            matrix = inMatrix;
        } else {
            matrix = formXtXInverse(inMatrix);
        }

        try {
            cholesky = (new CholeskyDecomposition(matrix)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in mvnOperator");
        }
    }

    public AdaptiveMetropolisOperator(Parameter parameter, double scaleFactor,
                                      MatrixParameter varMatrix, double weight, double beta, CoercionMode mode, boolean isVariance) {
        this(parameter, scaleFactor, varMatrix.getParameterAsMatrix(), weight, beta, mode, isVariance);
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
    	
    private double calculateCovariance (int number, double currentMatrixEntry, double[] values, int firstIndex, int secondIndex) {
    	
    	//number will always be > 1 here
    	double result  = currentMatrixEntry*(number - 1);
    	result += (values[firstIndex]*values[secondIndex]);
    	result += ((number-1)*oldMeans[firstIndex]*oldMeans[secondIndex] - number*newMeans[firstIndex]*newMeans[secondIndex]);
    	result /= ((double)number);
    	
    	return result;
    	
    }

    public double doOperation() throws OperatorFailedException {
    	
    	iterations++;
    	System.err.println("Using adaptive Metropolis (AM) operator: " + iterations);
    	
    	double[] x = parameter.getParameterValues();
    	
		double[] oldX = new double[x.length];
        System.arraycopy(x, 0, oldX, 0, x.length);
    	
    	if (iterations > 1) {
    		
    		//first recalculate the means using recursion
    		for (int i = 0; i < dim; i++) {
    			newMeans[i] = ((oldMeans[i]*(iterations-1)) + x[i])/iterations;
    		}
    		
    		//here we can simply use the double[][] matrix
    		for (int i = 0; i < dim; i++) {
    			for (int j = 0; j < dim; j++) {
    				empirical[i][j] = calculateCovariance(iterations, empirical[i][j], x, i, j);
    			}
    		}
    		
    		//test routine: first update old and new means
    		/*if (iterations == 2) {
    			System.err.println("old means:");
    			for (int i = 0; i < 5; i++) {
    				System.err.println(oldMeans[i]);
    			}
    			System.err.println("new means:");
    			for (int i = 0; i < 5; i++) {
    				System.err.println(newMeans[i]);
    			}
    			System.err.println("new values:");
    			for (int i = 0; i < 5; i++) {
    				System.err.println(x[i]);
    			}
    			System.err.println("empirical covariance matrix:");
    			for (int i = 0; i < 5; i++) {
    				for (int j = 0; j < 5; j++) {
    					System.err.print(empirical[i][j] + " ");
    				}
    				System.err.println();
    			}
    		}
    		System.exit(0);*/
    		
    	} else {
    		
    		//iterations == 1
    		for (int i = 0; i < dim; i++) {
        		oldMeans[i] = x[i];
        		newMeans[i] = x[i];
        	}
    		
    		for (int i = 0; i < dim; i++) {
    			for (int j = 0; j < dim; j++) {
    				empirical[i][j] = 0.0;
    			}
    		}
    	}
    	
    	if (iterations > cutoff) {
    		
    		System.err.println("Using empirical covariance matrix");
    		
    		/*double[] oldX = new double[x.length];
            System.arraycopy(x, 0, oldX, 0, x.length);*/
            
            double[] epsilon = new double[dim];
            
            for (int i = 0; i < dim; i++) {
                epsilon[i] = scaleFactor * MathUtils.nextGaussian();
            }
            
            double[][] proposal = new double[dim][dim];
            for (int i = 0; i < dim; i++) {
            	for (int j = 0; j < dim; j++) {
            		proposal[i][j] = (1 - beta)*Math.pow(2.38, 2)*empirical[i][j]/((double)dim) + beta*matrix[i][j];
            	}
            }
            
            //not necessary for first test phase, but will need to be performed when covariance matrix is being updated
        	try {
                cholesky = (new CholeskyDecomposition(proposal)).getL();
            } catch (IllegalDimension illegalDimension) {
                throw new RuntimeException("Unable to decompose matrix in MultivariateNormalMixtureOperator");
            }
            
            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    x[i] += cholesky[j][i] * epsilon[j];
                    // caution: decomposition returns lower triangular
                }
                parameter.setParameterValue(i, x[i]);
            }
            
            /*for (int i = 0; i < dim; i++) {
            	System.err.println(oldX[i] + " -> " + parameter.getValue(i));
            }*/
            
            //System.exit(0);
    		
    	} else {
    		
    		System.err.println("Using initial covariance matrix");
    		
            /*double[] oldX = new double[x.length];
            System.arraycopy(x, 0, oldX, 0, x.length);*/
            
            double[] epsilon = new double[dim];

            for (int i = 0; i < dim; i++) {
                epsilon[i] = scaleFactor * MathUtils.nextGaussian();
            }

            for (int i = 0; i < dim; i++) {
                for (int j = i; j < dim; j++) {
                    x[i] += cholesky[j][i] * epsilon[j];
                    // caution: decomposition returns lower triangular
                }
                parameter.setParameterValue(i, x[i]);
            }
            
            /*for (int i = 0; i < dim; i++) {
            	System.err.println(oldX[i] + " -> " + parameter.getValue(i));
            }*/
    		
    	}
    	
    	//copy new means to old means for next update iteration
    	System.arraycopy(newMeans, 0, oldMeans, 0, dim);
    	
        return 0;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "adapativeMetropolis(" + parameter.getParameterName() + ")";
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
            return AM_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        	//System.err.println("Parsing adaptive Metropolis (AM) operator.");
        	
            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double beta = xo.getDoubleAttribute(BETA);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

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
            	init[i].setParameterValue(i, Math.pow(0.1, 2)/((double)dim));
            }
            MatrixParameter varMatrix = new MatrixParameter(null, init);
            
            /*double[][] test = varMatrix.getParameterAsMatrix();
            for (int i = 0; i < dim; i++) {
            	for (int j = 0; j < dim; j++) {
            		System.err.print(test[i][j] + " ");
            	}
            	System.err.println();
            }*/
            
            // Make sure varMatrix is square and dim(varMatrix) = dim(parameter)

            if (!formXtXInverse) {
                if (varMatrix.getColumnDimension() != varMatrix.getRowDimension())
                    throw new XMLParseException("The variance matrix is not square");
            }

            if (varMatrix.getColumnDimension() != parameter.getDimension())
                throw new XMLParseException("The parameter and variance matrix have differing dimensions");

            return new AdaptiveMetropolisOperator(parameter, scaleFactor, varMatrix, weight, beta, mode, !formXtXInverse);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an adaptive Metropolis (AM) operator on a given parameter.";
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
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(FORM_XTX, true),
                new ElementRule(Parameter.class)
        };

    };
}
