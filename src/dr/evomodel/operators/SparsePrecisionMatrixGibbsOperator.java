/*
 * SparsePrecisionMatrixGibbsOperator.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.operators;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.GWishartLikelihood;//-------------------- Felipe
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.distribution.WishartGammalDistributionModel;
import dr.inference.model.DiagonalConstrainedMatrixView;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;//-------------------- Felipe
import dr.math.distributions.GWishartStatistics;//-------------------- Felipe
import dr.math.distributions.WishartDistribution;
import dr.math.distributions.WishartStatistics;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Attribute;
import dr.xml.*;

import java.util.List;
import java.util.Arrays; //-------------------- Felipe

/**
 * @author Marc Suchard
 * @author Felipe Pinheiro
 */
public class SparsePrecisionMatrixGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private static final String VARIANCE_OPERATOR = "sparsePrecisionGibbsOperator";
    public static final String TREE_MODEL = "treeModel";
    public static final String DISTRIBUTION = "distribution";
    public static final String PRIOR = "prior";
    private static final String WORKING = "workingDistribution";

    private final AbstractMultivariateTraitLikelihood traitModel;
    private AbstractMultivariateTraitLikelihood debugModel = null;
    private final ConjugateWishartStatisticsProvider conjugateWishartProvider;
    private final MultivariateDistributionLikelihood multivariateLikelihood;
    private final Parameter meanParam;
    private final MatrixParameterInterface precisionParam;

    private Statistics priorStatistics;
    private Statistics workingStatistics;
    
	private final GWishartStatistics priorGWishart;//------------------------ Felipe
	private boolean gWishartIsPrior = false;//------------------------ Felipe
	private Parameter priorAdjParameter;//------------------------ Felipe
	private double[][]priorAdjMatrix;//------------------------ Felipe
    
	private double priorDf;
    private SymmetricMatrix priorInverseScaleMatrix;

    private final MutableTreeModel treeModel;
    private final int dim;
    private double numberObservations;
    private final String traitName;
    private final boolean isSampledTraitLikelihood;
    private double pathWeight = 1.0;

    private boolean wishartIsModel = false;
    private WishartGammalDistributionModel priorModel = null;

    public SparsePrecisionMatrixGibbsOperator(
            MultivariateDistributionLikelihood likelihood,
            WishartStatistics priorDistribution,
            double weight) {
        super();
System.out.println("-------------------------------------------------------------- SparsePrecisionMatrixGibbsOperator C1");
        // Unnecessary variables
        this.traitModel = null;
        this.treeModel = null;
        this.traitName = null;
        this.conjugateWishartProvider = null;
        this.isSampledTraitLikelihood = false;
		this.priorGWishart=null;//------------------------------------ Felipe

        this.multivariateLikelihood = likelihood;
        MultivariateNormalDistributionModel density = (MultivariateNormalDistributionModel) likelihood.getDistribution();
        this.meanParam = density.getMeanParameter();
        this.precisionParam = density.getPrecisionMatrixParameter();
        this.dim = meanParam.getDimension();

        setupWishartStatistics(priorDistribution); // TODO Deprecate
        priorStatistics = setupStatistics(priorDistribution);

        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        setWeight(weight);
    }

    private class Statistics {
        final double degreesOfFreedom;
        final double[][] rateMatrix;
		final double[][] adjMatrix; //------------------------------Felipe

        Statistics(double degreesOfFreedom, double[][] rateMatrix) {
            this.degreesOfFreedom = degreesOfFreedom;
            this.rateMatrix = rateMatrix;
			this.adjMatrix=null; //------------------------------Felipe
        }
        
        Statistics(double degreesOfFreedom, double[][] rateMatrix, double[][] adjMatrix) {//------------------------------Felipe
			this.degreesOfFreedom = degreesOfFreedom;//------------------------------Felipe
			this.rateMatrix = rateMatrix;//------------------------------Felipe
			this.adjMatrix=adjMatrix;//------------------------------Felipe
		}
    }

    private Statistics setupStatistics(WishartStatistics distribution) {

        double[][] scale = distribution.getScaleMatrix();
        double[][] rate = null;

        if (scale != null) {
            rate = (new SymmetricMatrix(scale)).inverse().toComponents();
        }

        return new Statistics(distribution.getDF(), rate);
    }
    
    //------------------------------------------------------------------------------------------------FELIPE
    private Statistics setupStatistics(GWishartStatistics distribution) {    

		double[][] scale = distribution.getScaleMatrix(); 
		double[][] rate = null;

		if (scale != null) {
			rate = (new SymmetricMatrix(scale)).inverse().toComponents();
		}

		return new Statistics(distribution.getDF(), rate, distribution.getAdjMatrix());
	}
    
     
    private void setupWishartStatistics(WishartStatistics priorDistribution) {
        this.priorDf = priorDistribution.getDF();
        this.priorInverseScaleMatrix = null;
        double[][] scale = priorDistribution.getScaleMatrix();
        if (scale != null)
            this.priorInverseScaleMatrix =
                    (SymmetricMatrix) (new SymmetricMatrix(scale)).inverse();
    }
    
    //------------------------------------------------------------------------------------------------FELIPE
    private void setupGWishartStatistics(GWishartStatistics priorDistribution) {
		this.priorDf = priorDistribution.getDF();
		this.priorInverseScaleMatrix = null;
		double[][] scale = priorDistribution.getScaleMatrix();
		if (scale != null)
			this.priorInverseScaleMatrix = (SymmetricMatrix) (new SymmetricMatrix(scale)).inverse();
		this.priorAdjParameter=priorDistribution.getAdjParameter();
		this.priorAdjMatrix=priorDistribution.getAdjMatrix();
	}
    
    
    
    @Deprecated
    public SparsePrecisionMatrixGibbsOperator(
            MatrixParameterInterface precisionParam,
            AbstractMultivariateTraitLikelihood traitModel,
            WishartStatistics priorDistribution,
            double weight) {
        super();
        System.out.println("-------------------------------------------------------------- SparsePrecisionMatrixGibbsOperator C2");

        this.traitModel = traitModel;
        this.conjugateWishartProvider = null;
        this.meanParam = null;
        this.precisionParam = precisionParam;
		this.priorGWishart=null;//-----------------------------------------Felipe

        setupWishartStatistics(priorDistribution); // TODO Deprecate
        priorStatistics = setupStatistics(priorDistribution);

        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        setWeight(weight);

        this.treeModel = traitModel.getTreeModel();
        traitName = traitModel.getTraitName();
        dim = precisionParam.getRowDimension(); // assumed to be square

        isSampledTraitLikelihood = (traitModel instanceof SampledMultivariateTraitLikelihood);

        if (!isSampledTraitLikelihood &&
                !(traitModel instanceof ConjugateWishartStatisticsProvider)) {
            throw new RuntimeException("Only implemented for a SampledMultivariateTraitLikelihood or " +
                    "ConjugateWishartStatisticsProvider");
        }

        multivariateLikelihood = null;
    }

    
    //------------------------------------------------------------------------------------------------FELIPE
	@Deprecated
	public SparsePrecisionMatrixGibbsOperator(
			MatrixParameterInterface precisionParam,
			AbstractMultivariateTraitLikelihood traitModel, 
			GWishartStatistics priorDistribution, 
			double weight) {
		super();
		System.out.println("-------------------------------------------------------------- SparsePrecisionMatrixGibbsOperator C3");

		this.traitModel = traitModel;
		this.conjugateWishartProvider = null;
		this.meanParam = null;
		this.precisionParam = precisionParam;

		this.priorGWishart=priorDistribution;
		gWishartIsPrior = true;
		
		setupGWishartStatistics(priorDistribution); // TODO Deprecate
		priorStatistics = setupStatistics(priorDistribution);

		if (priorDistribution instanceof WishartGammalDistributionModel) {
			wishartIsModel = true;
			priorModel = (WishartGammalDistributionModel) priorDistribution;
		}

		setWeight(weight);

		this.treeModel = traitModel.getTreeModel();
		traitName = traitModel.getTraitName();
		dim = precisionParam.getRowDimension(); // assumed to be square

		isSampledTraitLikelihood = (traitModel instanceof SampledMultivariateTraitLikelihood);
		System.out.println("trait models inst of SampledMultivariateTraitLikelihood="
				+ (traitModel instanceof SampledMultivariateTraitLikelihood));
		System.out.println("trait models inst of ConjugateWishartStatisticsProvider="
				+ (traitModel instanceof ConjugateWishartStatisticsProvider));

		if (!isSampledTraitLikelihood && !(traitModel instanceof ConjugateWishartStatisticsProvider)) {
			throw new RuntimeException("Only implemented for a SampledMultivariateTraitLikelihood or "
					+ "ConjugateWishartStatisticsProvider");
		}

		multivariateLikelihood = null;
	}
	
    public SparsePrecisionMatrixGibbsOperator(
            ConjugateWishartStatisticsProvider wishartStatisticsProvider,
            MatrixParameterInterface extraPrecisionParam,
            WishartStatistics priorDistribution,
            WishartStatistics workingDistribution,
            double weight,
            AbstractMultivariateTraitLikelihood debugModel) {
        super();
        System.out.println("-------------------------------------------------------------- SparsePrecisionMatrixGibbsOperator C4");

        this.traitModel = null;
        this.debugModel = debugModel;
        this.conjugateWishartProvider = wishartStatisticsProvider;
        this.meanParam = null;
        this.precisionParam = (extraPrecisionParam != null ? extraPrecisionParam :
                conjugateWishartProvider.getPrecisionParameter());
        isSampledTraitLikelihood = false;
        this.treeModel = null;
        this.traitName = null;
		this.priorGWishart=null;//-----------------------------------------------Felipe


        setupWishartStatistics(priorDistribution); // TODO Deprecate
        priorStatistics = setupStatistics(priorDistribution);

        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        if (workingDistribution != null) {
            workingStatistics = setupStatistics(workingDistribution);
        }

        setWeight(weight);

        dim = precisionParam.getRowDimension(); // assumed to be square


        multivariateLikelihood = null;
    }
    
    //------------------------------------------------------------------------------------------------FELIPE
    public SparsePrecisionMatrixGibbsOperator(
            ConjugateWishartStatisticsProvider wishartStatisticsProvider,
            MatrixParameterInterface extraPrecisionParam,
            GWishartStatistics priorDistribution,
            WishartStatistics workingDistribution,
            double weight,
            AbstractMultivariateTraitLikelihood debugModel) {
        super();
        System.out.println("-------------------------------------------------------------- SparsePrecisionMatrixGibbsOperator C5");

        this.traitModel = null;
        this.debugModel = debugModel;
        this.conjugateWishartProvider = wishartStatisticsProvider;
        this.meanParam = null;
        this.precisionParam = (extraPrecisionParam != null ? extraPrecisionParam :
                conjugateWishartProvider.getPrecisionParameter());
        isSampledTraitLikelihood = false;
        this.treeModel = null;
        this.traitName = null;

		this.priorGWishart=priorDistribution;
		gWishartIsPrior = true;
		
		setupGWishartStatistics(priorDistribution); // TODO Deprecate
		priorStatistics = setupStatistics(priorDistribution);


        if (priorDistribution instanceof WishartGammalDistributionModel) {
            wishartIsModel = true;
            priorModel = (WishartGammalDistributionModel) priorDistribution;
        }

        if (workingDistribution != null) {
            workingStatistics = setupStatistics(workingDistribution);
        }

        setWeight(weight);

        dim = precisionParam.getRowDimension(); // assumed to be square


        multivariateLikelihood = null;
    }

    public void setPathParameter(double beta) {
        if (beta < 0 || beta > 1) {
            throw new IllegalArgumentException("Illegal path weight of " + beta);
        }
        pathWeight = beta;
    }

    public int getStepCount() {
        return 1;
    }

//    private void incrementScaledSquareMatrix(double[][] out, double[][] in, double scalar, int dim) {
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//                out[i][j] += scalar * in[i][j];
//            }
//        }
//    }

//    private void zeroSquareMatrix(double[][] out, int dim) {
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//                out[i][j] = 0.0;
//            }
//        }
//    }

    private void incrementOuterProduct(double[][] S,
                                       MultivariateDistributionLikelihood likelihood) {
    	System.out.println("incrementOuterProduct MultivariateDistributionLikelihood");
        double[] mean = likelihood.getDistribution().getMean();

        numberObservations = 0;
        List<Attribute<double[]>> dataList = likelihood.getDataList();

        for (Attribute<double[]> d : dataList) {

            double[] data = d.getAttributeValue();
            for (int i = 0; i < dim; i++) {
                data[i] -= mean[i];
            }

            for (int i = 0; i < dim; i++) {  // symmetric matrix,
                for (int j = i; j < dim; j++) {
                    S[j][i] = S[i][j] += data[i] * data[j];
                }
            }
            numberObservations += 1;
        }
    }

    private void incrementOuterProduct(double[][] S,
                                       ConjugateWishartStatisticsProvider integratedLikelihood) {
    	System.out.println("incrementOuterProduct ConjugateWishartStatisticsProvider");


        final WishartSufficientStatistics sufficientStatistics = integratedLikelihood.getWishartStatistics();
        final double[] outerProducts = sufficientStatistics.getScaleMatrix();

        final double df = sufficientStatistics.getDf();

        if (DEBUG) {
            System.err.println("OP df = " + df);
            System.err.println("OP    = " + new Vector(outerProducts));
        }

        if (debugModel != null) {
            final WishartSufficientStatistics debug = ((ConjugateWishartStatisticsProvider) debugModel).getWishartStatistics();
            System.err.println(df + " ?= " + debug.getDf());
            System.err.println(new Vector(outerProducts));
            System.err.println("");
            System.err.println(new Vector(debug.getScaleMatrix()));
            System.exit(-1);
        }

        final int dim = S.length;
        for (int i = 0; i < dim; i++) {
            System.arraycopy(outerProducts, i * dim, S[i], 0, dim);
        }
        numberObservations = df;


//        checkDiagonals(outerProducts);


    }

//    private void checkDiagonals(double[][] S) {
//        for (int i = 0; i < S.length; ++i) {
//            if (S[i][i] < 0.0) {
//                System.err.println("ERROR diag(S)\n" + new Matrix(S));
//                System.exit(-1);
//            }
//        }
//    } 

    private void incrementOuterProduct(double[][] S, NodeRef node) {
    	System.out.println("incrementOuterProduct NodeRef");

        if (!treeModel.isRoot(node)) {

            NodeRef parent = treeModel.getParent(node);
            double[] parentTrait = treeModel.getMultivariateNodeTrait(parent, traitName);
            double[] childTrait = treeModel.getMultivariateNodeTrait(node, traitName);
            double time = traitModel.getRescaledBranchLengthForPrecision(node);

            if (time > 0) {

                double sqrtTime = Math.sqrt(time);

                double[] delta = new double[dim];

                for (int i = 0; i < dim; i++)
                    delta[i] = (childTrait[i] - parentTrait[i]) / sqrtTime;

                for (int i = 0; i < dim; i++) {            // symmetric matrix,
                    for (int j = i; j < dim; j++)
                        S[j][i] = S[i][j] += delta[i] * delta[j];
                }
                numberObservations += 1; // This assumes a *single* observation per tip
            }
        }
        // recurse down tree
        for (int i = 0; i < treeModel.getChildCount(node); i++)
            incrementOuterProduct(S, treeModel.getChild(node, i));
    }

    private double[][] getOperationScaleMatrixAndSetObservationCount() {

        // calculate sum-of-the-weighted-squares matrix over tree
        double[][] S = new double[dim][dim];
        SymmetricMatrix S2;
        SymmetricMatrix inverseS2 = null;
        numberObservations = 0; // Need to reset, as incrementOuterProduct can be recursive

        if (isSampledTraitLikelihood) {
            incrementOuterProduct(S, treeModel.getRoot());
        } else { // IntegratedTraitLikelihood
            if (traitModel != null) { // is a tree
                incrementOuterProduct(S, (ConjugateWishartStatisticsProvider) traitModel); // TODO deprecate usage
            } else if (conjugateWishartProvider != null) {
                incrementOuterProduct(S, conjugateWishartProvider);
            } else { // is a normal-normal-wishart model
                incrementOuterProduct(S, multivariateLikelihood);
            }
        }

        try {
            S2 = new SymmetricMatrix(S);
            if (pathWeight != 1.0) {
                S2 = (SymmetricMatrix) S2.product(pathWeight);
            }
            if (priorInverseScaleMatrix != null)
                S2 = priorInverseScaleMatrix.add(S2);
            inverseS2 = (SymmetricMatrix) S2.inverse();

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        assert inverseS2 != null;

        return inverseS2.toComponents();
    }

    public double doOperation() {

        doOperationDontFireChange();
        precisionParam.fireParameterChangedEvent();

        return 0;
    }

    
    
    
    //------------------------------------------------------------------------------------------------FELIPE
	public void flipOne(int pos, Parameter parameter) {
        parameter.setParameterValue(pos, 0.0);
    }

    //------------------------------------------------------------------------------------------------FELIPE
    public void flipZero(int pos, Parameter parameter) {
        parameter.setParameterValue(pos, 1.0);
    }
    
    
    
        
    public void doOperationDontFireChange() {
        if (wishartIsModel) {
            setupWishartStatistics(priorModel); // TODO Deprecate
            priorStatistics = setupStatistics(priorModel);
        }

        final double[][] scaleMatrix = getOperationScaleMatrixAndSetObservationCount();
        final double treeDf = numberObservations;

        final double df = priorDf + treeDf * pathWeight;

        //------------------------------------------------------------------------------------------------FELIPE
        double[][] draw=null;
		
		if(gWishartIsPrior) { // draw will depend on the distribution if Wishart or G-Wishart ------------Felipe
				
			final int pos = MathUtils.nextInt(dim*(dim-1)/2);
			Parameter newAdjParameter= new Parameter.Default(priorAdjParameter.getParameterValues());

			int value = (int) newAdjParameter.getParameterValue(pos);
			
			System.out.println("\nGc0="+Arrays.toString(priorAdjParameter.getParameterValues()));
			System.out.println("Gp0="+Arrays.toString(newAdjParameter.getParameterValues()));

	        if (value == 0) {
	            flipZero(pos, newAdjParameter);
	        } else if (value == 1) {
	            flipOne(pos, newAdjParameter);
	        } else {
	            throw new RuntimeException("expected 1 or 0");
	        }
	        
	        System.out.println("Gp1="+Arrays.toString(newAdjParameter.getParameterValues()));
	        
	        //double[][] priorAdjMatrix = SymmetricMatrix.compoundSymmetricMatrix(new double[dim], priorAdjParameter.getParameterValues(), dim).toComponents();
	        //double[][] newAdjMatrix = SymmetricMatrix.compoundSymmetricMatrix(new double[dim], newAdjParameter.getParameterValues(), dim).toComponents();
	        
	        //System.out.println("priorAdjMatrix="+ Arrays.deepToString(priorAdjMatrix));
	        //System.out.println("  newAdjMatrix="+ Arrays.deepToString(newAdjMatrix));
	        
	        System.out.println("priorscaleMatrix="+ Arrays.deepToString(priorGWishart.getScaleMatrix()));
	        System.out.println("     scaleMatrix="+ Arrays.deepToString(scaleMatrix));
	        
	        double IGcPrior = GWishartLikelihood.computeNormalizationConstant(priorDf, priorGWishart.getScaleMatrix(), priorAdjParameter, dim);
			double IGcPosterior = GWishartLikelihood.computeNormalizationConstant(df, scaleMatrix, priorAdjParameter, dim);
			double IGpPrior = GWishartLikelihood.computeNormalizationConstant(priorDf, priorGWishart.getScaleMatrix(), newAdjParameter, dim);
			double IGpPosterior = GWishartLikelihood.computeNormalizationConstant(df, scaleMatrix, newAdjParameter, dim);

			System.out.println("IGcPrior="+IGcPrior);
			System.out.println("IGcPosterior="+IGcPosterior);
			System.out.println("IGpPrior="+IGpPrior);
			System.out.println("IGpPosterior="+IGpPosterior);
			  
			double alpha = 0.0;
			if (Double.isInfinite(IGcPrior) || Double.isInfinite(IGcPosterior) || Double.isInfinite(IGpPrior) || Double.isInfinite(IGpPosterior)) {
				alpha = 0.0;
	        }else {
	        	alpha = Math.exp(-(IGpPosterior - IGpPrior)  + (IGcPosterior - IGcPrior) );
	        }    
			

			double rand = MathUtils.nextDouble();			

		    if( rand<alpha && !Double.isNaN(alpha)) {
				System.out.println("U="+rand+" < alpha="+alpha+"--------------------------------------------------------------- ACCEPT NEW GRAPH");

		    	value = (int) priorAdjParameter.getParameterValue(pos);  
		    	if (value == 0) {
		    		flipZero(pos, priorAdjParameter);
			    } else if (value == 1) {
			        flipOne(pos, priorAdjParameter);
			    } else {
			        throw new RuntimeException("expected 1 or 0");
			    }
		    }else {
				System.out.println("U="+rand+" >= alpha="+alpha+"------------------------------------------- REJECT NEW GRAPH");

		    	value = (int) newAdjParameter.getParameterValue(pos);
		    	if (value == 0) {
			          flipZero(pos, newAdjParameter);
			    } else if (value == 1) {
			          flipOne(pos, newAdjParameter);
			    } else {
			    	throw new RuntimeException("expected 1 or 0");
			    	}
		    }    


		    System.out.println("Gc1="+Arrays.toString(priorAdjParameter.getParameterValues()));
	        System.out.println("Gp2="+Arrays.toString(newAdjParameter.getParameterValues()));

	
		    draw = GWishartLikelihood.nextGWishart(df, scaleMatrix, priorGWishart.getAdjMatrix()); 	        
    
		}else {
			draw = WishartDistribution.nextWishart(df, scaleMatrix);
		}
        if (DEBUG) {
            System.err.println("draw = " + new Matrix(draw));
        }

        for (int i = 0; i < dim; i++) {
            Parameter column = precisionParam.getParameter(i);
            for (int j = 0; j < dim; j++)
                column.setParameterValueQuietly(j, draw[j][i]);
        }
    }

    public MatrixParameterInterface getPrecisionParam() {
        return precisionParam;
    }

    public ConjugateWishartStatisticsProvider getConjugateWishartProvider() {
        return conjugateWishartProvider;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return VARIANCE_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return VARIANCE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(WEIGHT);

            AbstractMultivariateTraitLikelihood traitModel = (AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);
            ConjugateWishartStatisticsProvider ws = (ConjugateWishartStatisticsProvider) xo.getChild(ConjugateWishartStatisticsProvider.class);
            if (ws == traitModel) {
                ws = null;
            }
            
            

			GWishartLikelihood prior = null; //--------------------Felipe
            MatrixParameterInterface precMatrix = null;
            MultivariateDistributionLikelihood likelihood = null;

            if (traitModel != null) {
                precMatrix = traitModel.getDiffusionModel().getPrecisionParameter();
				prior = (GWishartLikelihood) xo.getChild(GWishartLikelihood.class);//--------------------Felipe
            }

            if (ws != null) {
                precMatrix = ws.getPrecisionParameter();
				prior = (GWishartLikelihood) xo.getChild(GWishartLikelihood.class);//-------------Felipe
            }

            /*
            if (traitModel == null && ws == null) { // generic likelihood and prior
                for (int i = 0; i < xo.getChildCount(); ++i) {
                    MultivariateDistributionLikelihood density = (MultivariateDistributionLikelihood) xo.getChild(i);
                    if (density.getDistribution() instanceof WishartStatistics) {
                        prior = density;
                    } else if (density.getDistribution() instanceof MultivariateNormalDistributionModel) {
                        likelihood = density;
                        precMatrix = ((MultivariateNormalDistributionModel)
                                density.getDistribution()).getPrecisionMatrixParameter();
                    }
                }

                if (prior == null || likelihood == null) {
                    throw new XMLParseException(
                            "Must provide a multivariate normal likelihood and Wishart prior in element '" +
                                    xo.getName() + "'\n"
                    );
                }
            }*/

            if (!(prior.getDistribution() instanceof WishartStatistics)) {
				throw new XMLParseException("Only a Wishart or GWishart distribution are conjugate for Gibbs sampling");//---Felipe
            }

            // Make sure precMatrix is square and dim(precMatrix) = dim(parameter)
            if (precMatrix.getColumnDimension() != precMatrix.getRowDimension()) {
                throw new XMLParseException("The variance matrix is not square or of wrong dimension");
            }

            if (traitModel != null && ws == null) {

                if (precMatrix instanceof DiagonalConstrainedMatrixView) {
                    precMatrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

                    if (precMatrix == null) {
                        throw new XMLParseException("Must provide unconstrained precision matrix");
                    }
                }
                
                //---------------------------------------------------------------FELIPE
                if (prior.getDistribution() instanceof GWishartStatistics) {
    				return new SparsePrecisionMatrixGibbsOperator(precMatrix, traitModel,
    						(GWishartStatistics) prior.getDistribution(), weight);
    			}else {		
    				return new SparsePrecisionMatrixGibbsOperator(precMatrix, traitModel,
    						(WishartStatistics) prior.getDistribution(), weight);
    				}
                
            } else if (ws != null) {

                if (precMatrix instanceof DiagonalConstrainedMatrixView) {
                    precMatrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

                    if (precMatrix == null) {
                        throw new XMLParseException("Must provide unconstrained precision matrix");
                    }
                } else {
                    precMatrix = null;
                }

                WishartStatistics workingDistribution = null;
                if (xo.hasChildNamed(WORKING)) {
                    workingDistribution = (WishartStatistics) xo.getElementFirstChild(WORKING);
                }

               
                
                if (prior.getDistribution() instanceof GWishartStatistics) {
    				return new SparsePrecisionMatrixGibbsOperator(
    						ws, precMatrix,	(GWishartStatistics) prior.getDistribution(), 
    						workingDistribution,
    						weight, traitModel);
    			}else {		
    				 return new SparsePrecisionMatrixGibbsOperator(
    	                        ws, precMatrix, (WishartStatistics) prior.getDistribution(),
    	                        workingDistribution,
    	                        weight, traitModel
    	                );
    				}

            } else {
                return new SparsePrecisionMatrixGibbsOperator(likelihood, (WishartStatistics) prior.getDistribution(), weight);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate normal random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
        		AttributeRule.newDoubleRule(WEIGHT),
				new ElementRule(AbstractMultivariateTraitLikelihood.class, true),
				new ElementRule(ConjugateWishartStatisticsProvider.class, true),
				new ElementRule(MultivariateDistributionLikelihood.class, true),
				new ElementRule(GWishartLikelihood.class, true),
				new ElementRule(MatrixParameterInterface.class, true),
        };
    };

    private static final boolean DEBUG = false;
}
