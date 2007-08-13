package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.GMRFSkylineLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;
import no.uib.cipr.matrix.*;

/**
 * @author Marc Suchard
 * @author Erik Bloomquist
 *         <p/>
 *         Created by IntelliJ IDEA.
 *         User: ebloomqu
 *         Date: Aug 8, 2007
 *         Time: 10:41:17 AM
 *         To change this template use File | Settings | File Templates.
 */
public class GMRFSkylineBlockUpdateOperator extends SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String BLOCK_UPDATE_OPERATOR = "gmrfBlockUpdateOperator";

    private double scaleFactor;
    private double lambdaScaleFactor;
    private int mode = CoercableMCMCOperator.DEFAULT;
    private int weight = 1;
    private int fieldLength;

    private Parameter popSizeParameter;
    private Parameter precisionParameter;
    private Parameter lambdaParameter;

    GMRFSkylineLikelihood gmrfField;

    private DenseVector one;
    private DenseMatrix I;
    
    public GMRFSkylineBlockUpdateOperator(GMRFSkylineLikelihood gmrfLikelihood, int weight, int mode) {
        super();
        this.mode = mode;
        this.weight = weight;
        gmrfField = gmrfLikelihood;
        popSizeParameter = gmrfLikelihood.getPopSizeParameter();
        precisionParameter = gmrfLikelihood.getPrecisionParameter();
        lambdaParameter = gmrfLikelihood.getLambdaParameter();

        scaleFactor = 4;
        lambdaScaleFactor = 0.0;
        fieldLength = popSizeParameter.getDimension();


        one = new DenseVector(fieldLength);
        I = new DenseMatrix(fieldLength, fieldLength);

        for (int i = 0; i < fieldLength; i++) {
            one.set(i, 1);
            I.set(i, i, 1.0);
        }

    }

    private double getNewLambda(double currentValue, double lambdaScale){
		double a = MathUtils.nextDouble()*lambdaScale - lambdaScale/2;
		double b = currentValue + a;
		if(b > 1)
			b = 2 - b;
		if(b < 0)
			b = -b;
		
		return b;
	}
    
    private double getNewPrecision(double currentValue, double scaleFactor) {
        double length = scaleFactor - 1 / scaleFactor;
        double returnValue = currentValue;


        if (scaleFactor == 1)
            return currentValue;
        if (MathUtils.nextDouble() < length / (length + 2 * Math.log(scaleFactor))) {
            returnValue = (1 / scaleFactor + length * MathUtils.nextDouble()) * currentValue;
        } else {
            returnValue = Math.pow(scaleFactor, 2.0 * MathUtils.nextDouble() - 1) * currentValue;
        }

        return returnValue;
    }


    //public double precisionPrior(double precision, double priorA, double priorB) {
    //    return (priorA - 1.0) * Math.log(precision) - precision * priorB;
    //}


    public DenseVector getMultiNormal(DenseVector Mean, UpperSPDDenseMatrix Variance) {
        DenseVector tempValue = new DenseVector(fieldLength);
        DenseVector returnValue = new DenseVector(fieldLength);
        UpperSPDDenseMatrix ab = Variance.copy();

        for (int i = 0; i < returnValue.size(); i++)
            tempValue.set(i, MathUtils.nextGaussian());

        DenseCholesky chol = new DenseCholesky(fieldLength, true);
        chol.factor(ab);

        UpperTriangDenseMatrix x = chol.getU();

        x.transMult(tempValue, returnValue);
        returnValue.add(Mean);
        return returnValue;
    }

    public double logGeneralizedDeterminant(SymmTridiagMatrix X) {
        //Set up the eigenvalue solver
        SymmTridiagEVD eigen = new SymmTridiagEVD(X.numRows(), false);
        //Solve for the eigenvalues
        try {
            eigen.factor(X);
        } catch (NotConvergedException e) {
            throw new RuntimeException("Not converged error in generalized determinate calculation.\n" + e.getMessage());
        }

        //Get the eigenvalues
        double[] x = eigen.getEigenvalues();

        double a = 0;
        for (int i = 0; i < x.length; i++) {
            if (x[i] > 0.00000001)
                a += Math.log(x[i]);
        }
        return a;
    }

    public double logGeneralizedDeterminant(UpperSPDPackMatrix X) {
        //Set up the eigenvalue solver
        SymmPackEVD eigen = new SymmPackEVD(X.numRows(), false);
        //Solve for the eigenvalues
        try {
            eigen.factor(X);
        } catch (NotConvergedException e) {
            throw new RuntimeException("Not converged error in generalized determinate calculation.\n" + e.getMessage());
        }

        //Get the eigenvalues;
        double[] x = eigen.getEigenvalues();

        double a = 0;
        for (int i = 0; i < x.length; i++) {

            if (x[i] > 0.00000001)
                a += Math.log(x[i]);
        }
        return a;

    }

    public static DenseVector newtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ) {

        DenseVector iterateGamma = currentGamma.copy();
        int numberIterations = 0;
        int maxIterations = 200;
        //TODO Integrate maxIterations and the stopValue into BEAST 
        while (gradient(data, iterateGamma, proposedQ).norm(Vector.Norm.Two) > 0.01) {
            inverseJacobian(data, iterateGamma, proposedQ).multAdd(gradient(data, iterateGamma, proposedQ), iterateGamma);
            numberIterations++;
        }
        
        if(numberIterations > maxIterations)
        	throw new RuntimeException("Newton Raphson algorithm did not converge within " + maxIterations + " step.\n");
        
        return iterateGamma;
    }

    private static DenseVector gradient(double[] data, DenseVector value, SymmTridiagMatrix Q) {

        DenseVector returnValue = new DenseVector(value.size());
        Q.mult(value, returnValue);
        for (int i = 0; i < value.size(); i++) {
            returnValue.set(i, -returnValue.get(i) - 1 + data[i] * Math.exp(-value.get(i)));
        }
        return returnValue;
    }


    private static DenseMatrix inverseJacobian(double[] data, DenseVector value, SymmTridiagMatrix Q) {

        SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
        for (int i = 0; i < value.size(); i++) {
            jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data[i]);
        }

        DenseMatrix inverseJacobian = Matrices.identity(jacobian.numRows());
        jacobian.solve(Matrices.identity(value.size()), inverseJacobian);

        return inverseJacobian;
    }


    public double doOperation() throws OperatorFailedException {

        double currentPrecisionParameter = precisionParameter.getParameterValue(0);
        
        //Generate a new precision.
        double newPrecisionParameter =
                this.getNewPrecision(currentPrecisionParameter, scaleFactor);

        precisionParameter.setParameterValue(0, newPrecisionParameter);

        //Generate a new mixing parameter 
        double currentLambda = lambdaParameter.getParameterValue(0);
        double newLambda = getNewLambda(currentLambda,lambdaScaleFactor);
        
        lambdaParameter.setParameterValue(0, newLambda);
        
        //Conditional on tau, generate a new value for gamma (it takes some work to do this).
        DenseVector currentGamma = new DenseVector(gmrfField.getPopSizeParameter().getParameterValues());
        DenseVector proposedGamma;

        //Get the current Precision and New Precision
        SymmTridiagMatrix currentQ = gmrfField.getScaledWeightMatrix(currentPrecisionParameter, currentLambda);
        SymmTridiagMatrix proposedQ = gmrfField.getScaledWeightMatrix(newPrecisionParameter, newLambda);

        //Get the data
        double[] wNative = gmrfField.getSufficientStatistics();

        //All of these Matricies and vectors are involved in the proposal steps
        SymmTridiagMatrix forwardQW = new SymmTridiagMatrix(proposedQ, true);
        SymmTridiagMatrix backwardQW = new SymmTridiagMatrix(currentQ, true);

        DenseVector forwardMean = new DenseVector(fieldLength);
        DenseVector backwardMean = new DenseVector(fieldLength);

        //Just some temporary storage.
        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector diagonal2 = new DenseVector(fieldLength);
        DenseVector diagonal3 = new DenseVector(fieldLength);

        //Needed for the inverses
        DenseMatrix workingDenseMatrix = Matrices.identity(fieldLength);

        //This finds the mode under Q proposed
        //NOT NEEDED IF EXPANDING UNDER CURRENT LOCATION
        DenseVector modeForward = newtonRaphson(wNative, currentGamma, proposedQ.copy());

       
        //This part determines whether the taylor expansion
        //will occur around the current value or the mode
        for (int i = 0; i < fieldLength; i++) {
            //Option1: Expand around current location
            //diagonal1.set(i, w[i]*Math.exp(-currentGamma.get(i)));
            //diagonal2.set(i, currentGamma.get(i) + 1);

            //Option2: Expand about the mode
            diagonal1.set(i, wNative[i] * Math.exp(-modeForward.get(i)));
            diagonal2.set(i, modeForward.get(i) + 1);

            forwardQW.set(i, i, diagonal1.get(i) + forwardQW.get(i, i));
            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
        }

        //Find the mean and variance of proposal distribution
        forwardQW.solve(I, workingDenseMatrix);
        UpperSPDDenseMatrix inverseForwardQW = new UpperSPDDenseMatrix(workingDenseMatrix, true);

        inverseForwardQW.mult(diagonal1, forwardMean);

        //Once we obtain the mean and Variance of the proposal
        //distribution, get a proposal!
        proposedGamma = getMultiNormal(forwardMean, inverseForwardQW);
        
        for (int i = 0; i < fieldLength; i++)
            popSizeParameter.setParameterValue(i, proposedGamma.get(i));

        //Now do some Metropolis-Hastings stuff in 5 steps.
        double hRatio = 0;

     //     double proposedLike = 0;
     //     double currentLike = 0;

        //1. First calculate the difference in the log-likelihoods
     //   for (int i = 0; i < fieldLength; i++) {
     //       proposedLike += -proposedGamma.get(i) - wNative[i] * Math.exp(-proposedGamma.get(i));
     //       currentLike += -currentGamma.get(i) - wNative[i] * Math.exp(-currentGamma.get(i));
     //   }

     //   hRatio = proposedLike - currentLike;

     //   2. Now the prior for gamma
     //   diagonal1.zero();
     //   diagonal2.zero();

     //   proposedQ.mult(proposedGamma, diagonal1);
     //  currentQ.mult(currentGamma, diagonal2);

     //  hRatio += 0.5 * logGeneralizedDeterminant(proposedQ) - 0.5 * proposedGamma.dot(diagonal1);
     //  hRatio -= 0.5 * logGeneralizedDeterminant(currentQ) - 0.5 * currentGamma.dot(diagonal2);

    //    //3. Now the prior for tau, no need for lambda prior, its 1.
    //  hRatio += precisionPrior(newPrecisionParameter, 0.01, 0.01);
    //  hRatio -= precisionPrior(currentPrecisionParameter, 0.01, 0.01);

        //4. Next find the difference in the proposal ratio's for gamma
        //This part sucks cuz you have to do all the backwards stuff
        diagonal1.zero();
        diagonal2.zero();
        diagonal3.zero();
        workingDenseMatrix = Matrices.identity(fieldLength);

        //This finds the mode under Q proposed
        DenseVector modeBackward = newtonRaphson(wNative, proposedGamma, currentQ.copy());

        for (int i = 0; i < fieldLength; i++) {
            //Expand around current value
            //diagonal1.set(i,w[i]*Math.exp(-proposedGamma.get(i)));
            //diagonal2.set(i, proposedGamma.get(i) + 1);

            //Expand about the mode
            diagonal1.set(i, wNative[i] * Math.exp(-modeBackward.get(i)));
            diagonal2.set(i, modeBackward.get(i) + 1);

            backwardQW.set(i, i, diagonal1.get(i) + backwardQW.get(i, i));
            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
        }

        backwardQW.solve(I, workingDenseMatrix);
        UpperSPDDenseMatrix inverseBackwardQW = new UpperSPDDenseMatrix(workingDenseMatrix, true);


        inverseBackwardQW.mult(diagonal1, backwardMean);

        //A lot of these steps are necessary
        //because MTJ has some bad matrix multiplication
        //and addition routines.

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, currentGamma.get(i) - backwardMean.get(i));
            diagonal2.set(i, proposedGamma.get(i) - forwardMean.get(i));
        }

        backwardQW.mult(diagonal1, diagonal3);

        hRatio += 0.5 * logGeneralizedDeterminant(backwardQW) - 0.5 * diagonal1.dot(diagonal3);

        diagonal3.zero();
        forwardQW.mult(diagonal2, diagonal3);

        hRatio -= 0.5 * logGeneralizedDeterminant(forwardQW) - 0.5 * diagonal2.dot(diagonal3);

       


       
        return hRatio;
    }

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
