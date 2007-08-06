package dr.evomodel.coalescent.operators;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import dr.evomodel.coalescent.GMRFSkylineLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Aug 6, 2007
 * Time: 10:41:17 AM
 * To change this template use File | Settings | File Templates.
 */
public class GMRFSkylineBlockUpdateOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String BLOCK_UPDATE_OPERATOR = "gmrfBlockUpdateOperator";

    private double scaleFactor = 0.5;
    private int mode = CoercableMCMCOperator.DEFAULT;
    private int weight = 1;
    private int fieldLength;
    GMRFSkylineLikelihood gmrfField;
    Parameter popSizeParameter;
    Parameter precisionParameter;
    Parameter lambdaParameter;
    DoubleMatrix2D weightMatrix;

    public GMRFSkylineBlockUpdateOperator(GMRFSkylineLikelihood gmrfLikelihood, int weight, int mode) {
        super();
        this.mode = mode;
        this.weight = weight;
        gmrfField = gmrfLikelihood;
        popSizeParameter = gmrfLikelihood.getPopSizeParameter();
        precisionParameter = gmrfLikelihood.getPrecisionParameter();
        lambdaParameter = gmrfLikelihood.getLambdaParameter();
        weightMatrix = gmrfLikelihood.getWeightMatrix();
        scaleFactor = 2.0;
        fieldLength = popSizeParameter.getDimension();


        I = new DenseDoubleMatrix2D(fieldLength, fieldLength);
        one = new DenseDoubleMatrix1D(fieldLength);

        for (int i = 0; i < fieldLength; i++) {
            I.set(i, i, 1);
            one.set(i, 1);
        }

        solver = new Algebra();
    }


    public DoubleMatrix2D getUpdateQ(double newTau, double newLambda) {
        DoubleMatrix2D a = new DenseDoubleMatrix2D(fieldLength, fieldLength);

        for (int i = 0; i < fieldLength; i++) {
            a.set(i, i, newTau * (weightMatrix.get(i, i) * newLambda + 1 - newLambda));
        }

        for (int i = 0; i < fieldLength - 1; i++) {
            a.set(i, i + 1, weightMatrix.get(i, i + 1) * newTau * newLambda);
            a.set(i + 1, i, weightMatrix.get(i + 1, i) * newTau * newLambda);
        }
        a.set(0, 0, weightMatrix.get(0, 0) * newTau);
        a.set(fieldLength - 1, fieldLength - 1,
                weightMatrix.get(fieldLength - 1, fieldLength - 1) * newTau);

        return a;

    }

    public double doOperation() throws OperatorFailedException {

//   	public void update(){

        //First generate a new value for tau
        double tau = precisionParameter.getParameterValue(0);
        double lambda = lambdaParameter.getParameterValue(0);
        double newTau = this.getNewTau(tau, scaleFactor);
//		double newLambda = this.getNewLambda(lambdaParameter.getParameterValue(0),0.0);
        double newLambda = 0.999999;

//        System.err.println("Old tau = "+tau+" New tau = "+newTau);
//        System.err.println("Old lambda = "+lambda+" New lambda = "+lambda);

        precisionParameter.setParameterValue(0, newTau);
        lambdaParameter.setParameterValue(0, newLambda);

        //Now conditional on tau.new, propose a new value for gamma
        DoubleMatrix1D mean = new DenseDoubleMatrix1D(fieldLength);
        //parm.getMu();
        double[] wNative = gmrfField.getSufficientStatistics();
        DoubleMatrix2D w = new DenseDoubleMatrix2D(fieldLength, fieldLength);
        for (int i = 0; i < fieldLength; i++)
            w.set(i, i, wNative[i]);
        //raw.getW();
//        System.err.println("Weight matrix:"+w);

        DoubleMatrix1D currentGamma = new DenseDoubleMatrix1D(popSizeParameter.getParameterValues());
//                parm.getGammaVector();

//        System.err.println("Current gamma = "+currentGamma);

        DoubleMatrix1D proposedGamma = new DenseDoubleMatrix1D(fieldLength);


        DoubleMatrix2D qProposed = getUpdateQ(newTau, newLambda);

//        System.err.println("qProposed = "+qProposed);


        DoubleMatrix2D qCurrent = getUpdateQ(tau, newLambda);

//        System.err.println("qCurrent = "+qCurrent);

        DoubleMatrix2D qwInverseProposed;
        DoubleMatrix2D qwInverseBackward;

        DoubleMatrix2D temp1 = new DenseDoubleMatrix2D(fieldLength, fieldLength);
        DoubleMatrix2D temp2 = new DenseDoubleMatrix2D(fieldLength, fieldLength);
        DoubleMatrix2D temp3 = new DenseDoubleMatrix2D(fieldLength, fieldLength);
        DoubleMatrix1D temp4 = new DenseDoubleMatrix1D(fieldLength);
        DoubleMatrix1D temp5 = new DenseDoubleMatrix1D(fieldLength);

        for (int i = 0; i < fieldLength; i++) {
            temp1.set(i, i, w.get(i, i) * Math.exp(-currentGamma.get(i)));
            temp2.set(i, i, currentGamma.get(i) + 1);

        }

        qwInverseProposed = solver.inverse(add(qProposed, temp1, true));

        temp3 = add(solver.mult(temp1, temp2), I, false);
        temp4 = solver.mult(temp3, one);

        temp4 = addVector(solver.mult(qProposed, mean), temp4, true);

        temp4 = solver.mult(qwInverseProposed, temp4);


        proposedGamma = getMultiNormal(temp4, qwInverseProposed);

//        System.err.println("Proposed gamma = "+proposedGamma);

        for (int i = 0; i < fieldLength; i++)
            popSizeParameter.setParameterValue(i, proposedGamma.get(i));

        //Now do some Metropolis-Hastings stuff in 6 steps.
        double hRatio = 0;

        //1. First calculate the difference in log-likelihoods
        for (int i = 0; i < fieldLength; i++) {
            hRatio += -proposedGamma.get(i) - w.get(i, 1) * Math.exp(-proposedGamma.get(i));
            hRatio -= -currentGamma.get(i) - w.get(i, i) * Math.exp(-currentGamma.get(i));
        }

        //2. Now the difference in proposal distributions.
        hRatio += gammaPrior(newTau, 0.1, 0.1);
        hRatio -= gammaPrior(tau, 0.1, 0.1);

        //3. Next find the difference in the proposal ratio's for gamma
        //This part sucks cuz you have to do all the backwards stuff
        for (int i = 0; i < fieldLength; i++) {
            temp1.set(i, i, w.get(i, i) * Math.exp(-proposedGamma.get(i)));
            temp2.set(i, i, proposedGamma.get(i) + 1);
        }

        qwInverseBackward = solver.inverse(add(qCurrent, w, true));
        temp3 = add(solver.mult(temp1, temp2), I, false);
        temp5 = solver.mult(temp3, one);
        temp5 = addVector(solver.mult(qCurrent, mean), temp4, true);
        temp5 = solver.mult(qwInverseBackward, temp4);

        hRatio += -1 / 2 * solver.mult(addVector(currentGamma, temp5, false),
                solver.mult(qwInverseBackward, addVector(currentGamma, temp5, false)));
        hRatio -= -1 / 2 * solver.mult(addVector(proposedGamma, temp4, false),
                solver.mult(qwInverseProposed, addVector(proposedGamma, temp4, false)));

        //4. Don't forget to inlucde the determinants of Q and Q^\star
//        System.err.println("hRatioB = "+hRatio);
//        System.err.println("detC = "+solver.det(qCurrent));
        double logDetRatio = Math.log(solver.det(qCurrent)) - Math.log(solver.det(qProposed));
//        System.err.println("detRatio = "+logDetRatio);

        hRatio += logDetRatio;

        //5. We don't have to worry about the proposal distribution for tau
        //because it is symmetric.

//		 //6. Finally, run the usual HM acceptance step.
//		 if(rand.getUniform(0,1) < Math.exp(Math.min(0,hRatio))){
//        	 mon.blockUpdateSuccess(true);
//			 parm.updateQ(newTau,newLambda);
//        	 parm.updateGamma(proposedGamma);
//         }
//		 else{
//			 mon.blockUpdateSuccess(false);
//		 }

        //Done!
//	}
        System.err.println("hRatio = " + hRatio);

        return hRatio;
    }


    public DoubleMatrix1D getMultiNormal(DoubleMatrix1D mean, DoubleMatrix2D variance) {
        DoubleMatrix1D Z = new DenseDoubleMatrix1D(mean.size());
        for (int i = 0; i < Z.size(); i++)
            Z.set(i, MathUtils.nextGaussian());
//                    norm.nextDouble());

        CholeskyDecomposition cholesky = new CholeskyDecomposition(variance);
//		Algebra solver = new Algebra();

        return addMatrix(mean, solver.mult(cholesky.getL(), Z));

    }

    private static DoubleMatrix1D addMatrix(DoubleMatrix1D A, DoubleMatrix1D B) {
        DoubleMatrix1D returnValue = new SparseDoubleMatrix1D(A.size());
        for (int i = 0; i < A.size(); i++)
            returnValue.set(i, A.get(i) + B.get(i));
        return returnValue;
    }


    private Parameter parm;
    private Random rand;
//	private RawData raw;
//	private Settings set;
//	private SamplerMonitor mon;

    private DoubleMatrix2D I;
    private DoubleMatrix1D one;
    private Algebra solver;

//	public Blockupdate(Parameters inP, Random inR,RawData data,Settings inSet,SamplerMonitor inM){
//		parm = inP;
//		rand = inR;
//		raw = data;
//		set = inSet;
//
//		I = new DenseDoubleMatrix2D(fieldLength, fieldLength);
//		one = new DenseDoubleMatrix1D(fieldLength);
//
//		for(int i = 0; i < fieldLength; i++){
//			I.set(i, i, 1);
//			one.set(i,1);
//		}
//
//		solver = new Algebra();
//		mon = inM;
//	}

    private double getNewLambda(double currentValue, double lambdaScale) {
        double a = MathUtils.nextDouble() * lambdaScale - lambdaScale / 2.0;
        double b = currentValue + a;
        if (b == 1)
            b = 0.999999;
        if (b == 0)
            b = 0.000001;
        if (b > 1)
            b = 2 - b;
        if (b < 0)
            b = -b;

        return b;
    }

    private double getNewTau(double currentValue, double scaleFactor) {
        double length = scaleFactor - 1 / scaleFactor;
        double returnValue = currentValue;

        if (scaleFactor == 1)
            return currentValue;
        if (//rand.getUniform(0, 1)
                MathUtils.nextDouble()
                        < length / (length + 2 * Math.log(scaleFactor))) {
            returnValue = (1 / scaleFactor + length * MathUtils.nextDouble())
//                    rand.getUniform(0, 1))
                    * currentValue;
        } else {
            returnValue = Math.pow(scaleFactor, 2.0 *//rand.getUniform(0, 1)
                    MathUtils.nextDouble()
                    - 1) * currentValue;
        }
        return returnValue;
    }

    public static DoubleMatrix2D add(DoubleMatrix2D A, DoubleMatrix2D B, boolean add) {
        DoubleMatrix2D returnValue = new DenseDoubleMatrix2D(A.rows(), A.columns());
        if (add) {
            for (int i = 0; i < A.rows(); i++)
                for (int j = 0; j < B.rows(); j++)
                    returnValue.set(i, j, A.get(i, j) + B.get(i, j));
        } else {
            for (int i = 0; i < A.rows(); i++)
                for (int j = 0; j < B.rows(); j++)
                    returnValue.set(i, j, A.get(i, j) - B.get(i, j));
        }
        return returnValue;
    }

    public static DoubleMatrix1D addVector(DoubleMatrix1D A, DoubleMatrix1D B, boolean add) {
        DoubleMatrix1D returnValue = new SparseDoubleMatrix1D(A.size());
        if (add) {
            for (int i = 0; i < A.size(); i++)
                returnValue.set(i, A.get(i) + B.get(i));
        } else {
            for (int i = 0; i < A.size(); i++)
                returnValue.set(i, A.get(i) - B.get(i));
        }
        return returnValue;
    }

    public double gammaPrior(double tau, double priorA, double priorB) {
        return (priorA - 1.0) * Math.log(tau) - tau * priorB;
    }

    /*
    * This is the method that does all the work.
    */

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return BLOCK_UPDATE_OPERATOR;
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
        throw new RuntimeException("Tried to automized block update");

    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public int getMode() {
        return mode;
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

    public int getWeight() {
        return weight;
    }

    public void setWeight(int w) {
        weight = w;
    }


    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return BLOCK_UPDATE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            int mode = CoercableMCMCOperator.DEFAULT;

            if (xo.hasAttribute(AUTO_OPTIMIZE)) {
                if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
                    mode = CoercableMCMCOperator.COERCION_ON;
                } else {
                    mode = CoercableMCMCOperator.COERCION_OFF;
                }
            }

            int weight = xo.getIntegerAttribute(WEIGHT);
//			double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
//
//			if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
//				throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
//			}

//            XMLObject cxo = (XMLObject) xo.getChild(POPULATION_PARAMETER);
//            Parameter populationSizeParameter = (Parameter) cxo.getChild(Parameter.class);

//            cxo = (XMLObject) xo.getChild(PRECISION_PARAMETER);
//            Parameter precisionParameter = (Parameter) cxo.getChild(Parameter.class);

//			XMLObject cxo = (XMLObject) xo.getChild(PRECISION_PRIOR);
//            DistributionLikelihood precisionPrior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);

            GMRFSkylineLikelihood gmrfLikelihood = (GMRFSkylineLikelihood) xo.getChild(GMRFSkylineLikelihood.class);

            return new GMRFSkylineBlockUpdateOperator(//populationSizeParameter, precisionParameter,
                    gmrfLikelihood, weight, mode);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a GMRF block-update operator for the joint distribution of the population sizes and precision parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//				AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newIntegerRule(WEIGHT),
//				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
//				new ElementRule(PRECISION_PRIOR, new XMLSyntaxRule[]{
//                new ElementRule(DistributionLikelihood.class),
//				}),
//                new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
//                        new ElementRule(Parameter.class)
//                }),
//                new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
//                        new ElementRule(Parameter.class)
//                }),
                new ElementRule(GMRFSkylineLikelihood.class)
        };

    };


}
