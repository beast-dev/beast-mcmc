/*
 * GaussianProcessSkytrackLikelihood.java
 *
 * Copyright (c) 2002-2011 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent;

//import com.lowagie.text.Paragraph;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
//import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GaussianProcessSkytrackLikelihoodParser;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import no.uib.cipr.matrix.*;

import java.util.ArrayList;
import java.util.List;



/**
 * @author Vladimir Minin
 * @author Marc Suchard
 * @author Julia Palacios
 */

//For implementation,
public class GaussianProcessSkytrackLikelihood extends OldAbstractCoalescentLikelihood {

//    protected Parameter groupSizeParameter;
    public static final double LOG_TWO_TIMES_PI = 1.837877;
	protected Parameter precisionParameter;
    protected Parameter lambda_boundParameter;
    protected Parameter numGridPoints;
    protected Parameter lambdaParameter;    //prior for lambda_bound, will be used in operators only
    protected Parameter betaParameter;
    protected Parameter alphaParameter;
    protected Parameter GPtype;
    protected Parameter GPcounts;
    protected Parameter coalfactor;
    protected Parameter popSizeParameter;     //before called GPvalues
    protected Parameter changePoints;
    protected Parameter popValue;
    protected Parameter CoalCounts;
    protected Parameter numPoints;

//    protected double [] GPchangePoints;
//    protected double [] storedGPchangePoints;
    protected double [] GPcoalfactor;
//    protected double [] storedGPcoalfactor;
    protected double [] GPCoalInterval;

//    protected double [] storedcoalfactor;
//    protected int [] GPcounts;   //It changes values, no need to storage
//    protected int [] storedGPcounts;
    protected int [] CoalPosIndicator;
    protected double [] CoalTime;
    protected int numintervals;
    protected int numcoalpoints;
    protected double constlik;

//    Those that change size, they are initialized per tree, no need to store them
//    use as Parameter since they will be changing by operators
//    protected Parameter GPtimepoints;  //tree + latent
//    protected double GPintervalkey;         // membership that links with those that do not change in size
//    protected Parameter GPcoalfactor2;        // choose(k,2) depending on membership
      // 1 if observed, -1 if latent
//    protected int[] storedGPtype;

//    public double[] GPvalues;     //may need to change type: Parameter? didn't know how to work with it

    protected double logGPLikelihood;
//    protected double storedLogGPLikelihood;
    protected SymmTridiagMatrix weightMatrix;       //this now changes in dimension, no need to storage
//	protected MatrixParameter dMatrix;
	protected boolean rescaleByRootHeight;


    private static List<Tree> wrapTree(Tree tree) {
        List<Tree> treeList = new ArrayList<Tree>();
        treeList.add(tree);
        return treeList;
    }

    public GaussianProcessSkytrackLikelihood(Tree tree,
                                             Parameter precParameter,
                                             boolean rescaleByRootHeight, Parameter numGridPoints,  Parameter lambda_bound, Parameter lambda_parameter, Parameter popParameter, Parameter alpha_parameter, Parameter beta_parameter, Parameter change_points, Parameter GPtype, Parameter GPcounts, Parameter coalfactor, Parameter popValues, Parameter CoalCounts, Parameter numPoints) {
        this(wrapTree(tree),  precParameter, rescaleByRootHeight, numGridPoints, lambda_bound, lambda_parameter, popParameter, alpha_parameter, beta_parameter, change_points,GPtype,GPcounts,coalfactor,popValues,CoalCounts, numPoints);

    }


    public GaussianProcessSkytrackLikelihood(String name) {
		super(name);
	}




    public GaussianProcessSkytrackLikelihood(List<Tree> treeList,
                                             Parameter precParameter,
                                              boolean rescaleByRootHeight, Parameter numGridPoints, Parameter lambda_bound, Parameter lambda_parameter, Parameter popParameter, Parameter alpha_parameter, Parameter beta_parameter, Parameter change_points, Parameter GPtype, Parameter GPcounts, Parameter coalfactor, Parameter popValues, Parameter CoalCounts, Parameter numPoints) {
        super(GaussianProcessSkytrackLikelihoodParser.SKYTRACK_LIKELIHOOD);



                this.popSizeParameter = popParameter;
                this.popValue=popValues;
                this.changePoints=change_points;
                this.numPoints=numPoints;
//                this.groupSizeParameter = groupParameter;
                this.precisionParameter = precParameter;
                this.lambdaParameter = lambda_parameter;
                this.betaParameter = beta_parameter;
                this.alphaParameter=alpha_parameter;
//                this.dMatrix = dMatrix;
                this.rescaleByRootHeight = rescaleByRootHeight;
                this.numGridPoints = numGridPoints;
                this.lambda_boundParameter= lambda_bound;
                this.GPcounts=GPcounts;
                this.GPtype=GPtype;
                this.coalfactor=coalfactor;
                this.CoalCounts=CoalCounts;


//                addVariable(GPvalues);
                addVariable(precisionParameter);
                addVariable(popSizeParameter);
                addVariable(changePoints);
                addVariable(numPoints);
//                addVariable(popValue);
                addVariable(GPcounts);
//                addVariable(GPcoalfactor);
                addVariable(GPtype);

                addVariable(coalfactor);
                addVariable(lambda_boundParameter);
                addVariable(CoalCounts);



//                addVariable(lambdaParameter);
//                addVariable(lambda_boundParameter);
//                if (betaParameter != null) {
//                    addVariable(betaParameter);
//                }

                setTree(treeList);

                wrapSetupIntervals();

//        intervalCount = the size for constant vectors



//                int fieldLength = getCorrectFieldLength();
                numintervals= getIntervalCount();
                numcoalpoints=getCorrectFieldLength();

        GPcoalfactor = new double[numintervals];
        GPCoalInterval=new double[numcoalpoints];
        CoalPosIndicator= new int[numcoalpoints];
        CoalTime=new double[numcoalpoints];


//        storedGPcoalfactor = new double[numintervals];
        GPcounts.setDimension(numintervals);
        CoalCounts.setDimension(numcoalpoints);


//        storedGPcounts= new int[numintervals];
        GPtype.setDimension(numcoalpoints);
        numPoints.setParameterValue(0,numcoalpoints);
//        storedGPtype = new int[numcoalpoints];
        popSizeParameter.setDimension(numcoalpoints);

        int gridpoint= (int) numGridPoints.getParameterValue(0);
        popValue.setDimension(gridpoint);

        changePoints.setDimension(numcoalpoints);
        coalfactor.setDimension(numcoalpoints);
//        storedcoalfactor= new double[numcoalpoints];




                initializationReport();
                setupSufficientStatistics();


        setupGPvalues();



//              System.err.println(getLogLikelihood());


         }


// Methods that override existent methods


    protected void setTree(List<Tree> treeList) {
        if (treeList.size() != 1) {
             throw new RuntimeException("GP-based method only implemented for one tree");
        }
        this.tree = treeList.get(0);
        this.treesSet = null;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }


    protected void wrapSetupIntervals() {
        setupIntervals();
    }
//
//    protected int getCorrectFieldLength() {
//            return tree.getExternalNodeCount() - 1;
//    }


//
//


    //This is actually the Augmented loglikelihood for fixed genealogy. For sequence data directly
//    this becomes the coalescent point process prior on an augmented "tree"
    public double calculateLogLikelihood(Parameter Gfunction, Parameter latentCounts, Parameter eventType, Parameter upper_Bound, double [] Gfactor) {
        double upperBound = upper_Bound.getParameterValue(0);
//        System.err.println("Likelihood with "+getPopSizeParameter().getSize()+"and G-function"+eventType.getSize());

        logGPLikelihood=-upperBound*getConstlik();

        for (int i=0; i<latentCounts.getSize(); i++){
            if (Gfactor[i]>0) {
          logGPLikelihood+=latentCounts.getParameterValue(i)*Math.log(upperBound*Gfactor[i]);
            }
        }
        double[] currentGfunction = Gfunction.getParameterValues();
        for (int i=0; i<Gfunction.getSize();i++){
        logGPLikelihood+= -Math.log(1+Math.exp(-eventType.getParameterValue(i)*currentGfunction[i]));
        }
        return logGPLikelihood;
    }

//    protected double calculateLogCoalescentLikelihood() {
//
//		if (!intervalsKnown) {
//			// intervalsKnown -> false when handleModelChanged event occurs in super.
//			wrapSetupIntervals();
//			setupGMRFWeights();
//            intervalsKnown = true;
//		}
//
//		// Matrix operations taken from block update sampler to calculate data likelihood and field prior
//
//		double currentLike = 0;
//        double[] currentGamma = popSizeParameter.getParameterValues();
//
//		for (int i = 0; i < fieldLength; i++) {
//			currentLike += -currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
//		}
//
//		return currentLike;// + LogNormalDistribution.logPdf(Math.exp(popSizeParameter.getParameterValue(coalescentIntervals.length - 1)), mu, sigma);
//	return 0.0;
//    }

//

    public double getConstlik(){
        return constlik;
    }


//For fixed genealogy this contains the Augmented likelihood, the GP prior and prior on a the upper bound
	public double getLogLikelihood() {
//            System.err.println("getlog used");
		if (!likelihoodKnown) {
			logLikelihood =
              calculateLogLikelihood(popSizeParameter,GPcounts,GPtype,lambda_boundParameter,GPcoalfactor)+calculateLogGP()
                        +getLogPriorLambda(lambdaParameter.getParameterValue(0),0.01,lambda_boundParameter.getParameterValue(0));
			likelihoodKnown = true;
		}
		return logLikelihood;
//        return 0.0;
	}

//Calculates prior on g function
    protected double calculateLogGP() {

           if (!intervalsKnown) {
//               System.err.println("intervalsknown");
//               intervalsKnown -> false when handleModelChanged event occurs in super.
               wrapSetupIntervals();
//               setupQmatrix(precisionParameter.getParameterValue(0));
               intervalsKnown = true;
           }


           setupQmatrix(precisionParameter.getParameterValue(0));

           double currentLike;
           DenseVector diagonal1 = new DenseVector(popSizeParameter.getSize());
           DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());

           SymmTridiagMatrix currentQ = weightMatrix;
//        System.err.println("Q.matrix sizes"+currentQ.numRows()+"and"+currentQ.numColumns());
           currentQ.mult(currentGamma, diagonal1);

           currentLike = -0.5 * logGeneralizedDeterminant(currentQ) - 0.5 * currentGamma.dot(diagonal1) - 0.5 * (popSizeParameter.getSize() - 1) * LOG_TWO_TIMES_PI;

        return currentLike;
       }

//    Calculates logprior on Upper Bound
    private double getLogPriorLambda(double lambdaMean, double epsilon, double lambdaValue){
    double res;
    if (lambdaValue < lambdaMean) {res=epsilon*(1/lambdaMean);}
    else {res=Math.log(1-epsilon)*(1/lambdaMean)*Math.exp(-(1/lambdaMean)*(lambdaValue-lambdaMean)); }
    return res;
}


    //log pseudo-determinant
       public static double logGeneralizedDeterminant(SymmTridiagMatrix X) {
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
           for (double d : x) {
               if (d > 0.00001)
                   a += Math.log(d);
           }

           return a;
       }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type){
		likelihoodKnown = false;
        // Parameters (precision and popsizes do not change intervals or GMRF Q matrix (I DON'T UNDERSTAND THIS)
	}

//
//protected void restoreState() {
//		super.restoreState();
//		System.arraycopy(storedcoalfactor, 0, coalfactor, 0, storedcoalfactor.length);
//        System.arraycopy(storedGPtype,0,GPtype,0,storedGPtype.length);
//        System.arraycopy(storedGPcoalfactor,0,GPcoalfactor,0,storedGPcoalfactor.length);
//        System.arraycopy(storedGPcounts,0,GPcounts,0,storedGPcounts.length);
////		weightMatrix = storedWeightMatrix;
//        logGPLikelihood = storedLogGPLikelihood;
//    }
//
//
//
//protected void storeState() {
//		super.storeState();
//	 	System.arraycopy(GPtype, 0, storedGPtype, 0, GPtype.length);
//        System.arraycopy(GPcoalfactor,0,storedGPcoalfactor,0,GPcoalfactor.length);
//		System.arraycopy(coalfactor, 0, storedcoalfactor, 0, coalfactor.length);
//        System.arraycopy(GPcounts, 0, storedGPcounts,0,GPcounts.length);
//
////		storedWeightMatrix = weightMatrix.copy();
//        storedLogGPLikelihood = logGPLikelihood;
//	}
//                I don't understand this
       public String toString() {
        return getId() + "(" + Double.toString(getLogLikelihood()) + ")";
    }
////    private final Parameter latentPoints;
//
//    private final Parameter lambda_bound;
//
//
//
    public void initializationReport() {

		System.out.println("Creating a GP based estimation of effective population trajectories:");
//		System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
		System.out.println("\tIf you publish results using this model, please reference: Minin, Palacios, Suchard (XXXX), AAA");

	}

    public static void checkTree(TreeModel treeModel) {

            // todo Should only be run if there exists a zero-length interval

//        TreeModel treeModel = (TreeModel) tree;
            for (int i = 0; i < treeModel.getInternalNodeCount(); i++) {
                NodeRef node = treeModel.getInternalNode(i);
                if (node != treeModel.getRoot()) {
                    double parentHeight = treeModel.getNodeHeight(treeModel.getParent(node));
                    double childHeight0 = treeModel.getNodeHeight(treeModel.getChild(node, 0));
                    double childHeight1 = treeModel.getNodeHeight(treeModel.getChild(node, 1));
                    double maxChild = childHeight0;
                    if (childHeight1 > maxChild)
                        maxChild = childHeight1;
                    double newHeight = maxChild + MathUtils.nextDouble() * (parentHeight - maxChild);
                    treeModel.setNodeHeight(node, newHeight);
                }
            }
            treeModel.pushTreeChangedEvent();

        }


     //Sufficient Statistics for GP - coal+sampling

    protected void setupSufficientStatistics() {


		double length = 0.0;
        double prevLength=0.0;
        int countcoal = 0;
        constlik= 0;
		for (int i = 0; i < getIntervalCount(); i++) {

			  length += getInterval(i);

//              if (GPcounts.getSize()<=i){
//                  GPcounts.addDimension(i,0.0);
//              }  else {
              GPcounts.setParameterValue(i,0.0);

//              }
              GPcoalfactor[i] =getLineageCount(i)*(getLineageCount(i)-1.0) / 2.0;
              constlik+=GPcoalfactor[i]*getInterval(i);

//            System.err.println("i: "+i+"val: "+length+" type: "+getIntervalType(i)+" lineages: "+getLineageCount(i));
              if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                    GPcounts.setParameterValue(i,1.0);
                    GPtype.setParameterValue(countcoal,1.0);

                    CoalPosIndicator[countcoal]=i;
                    changePoints.setParameterValue(countcoal,length);
                    CoalCounts.setParameterValue(countcoal,0.0);
                    CoalTime[countcoal]=length;

                    GPCoalInterval[countcoal]=length-prevLength;

                    coalfactor.setParameterValue(countcoal,getLineageCount(i)*(getLineageCount(i)-1)/2.0);

                    countcoal++;
                    prevLength=length;
			        }


        }


    }


    protected int getCorrectFieldLength() {
        return tree.getExternalNodeCount() - 1;
    }

    protected void setupQmatrix(double precision) {
//                   System.err.println("changepoints size"+changePoints.getSize());

            //Set up the weight Matrix
            double trick=0.000001;
            double[] offdiag = new double[changePoints.getSize() - 1];
            double[] diag = new double[changePoints.getSize()];


             for (int i = 0; i < changePoints.getSize() - 1; i++) {
                    offdiag[i] = precision*(-1.0 / (changePoints.getParameterValue(i+1)-changePoints.getParameterValue(i)));
                 if (i<getCorrectFieldLength()-2){
                    diag[i+1]= -offdiag[i]+precision*(1.0/(changePoints.getParameterValue(i+2)-changePoints.getParameterValue(i+1))+trick);

                 }
                }
//              Diffuse prior correction - intrinsic
             //Take care of the endpoints
            diag[0] = -offdiag[0]+precision*trick;

            diag[getCorrectFieldLength() - 1] = -offdiag[getCorrectFieldLength() - 2]+precision*(trick);
            weightMatrix = new SymmTridiagMatrix(diag, offdiag);
        }



    protected void setupGPvalues() {

        setupQmatrix(precisionParameter.getParameterValue(0));
        int length = getCorrectFieldLength();
        DenseVector StandNorm = new DenseVector(length);
        DenseVector MultiNorm = new DenseVector(length);
        for (int i=0; i<length;i++){
            StandNorm.set(i,MathUtils.nextGaussian());
//            StandNorm.set(i,0.1);
                      }
        UpperSPDBandMatrix Qcurrent = new UpperSPDBandMatrix(weightMatrix, 1);
        BandCholesky U = new BandCholesky(length,1,true);
        U.factor(Qcurrent);
        UpperTriangBandMatrix CholeskyUpper = U.getU();

        CholeskyUpper.solve(StandNorm,MultiNorm);
        for (int i=0; i<length;i++){
//            popSizeParameter.setParameterValue(i,MultiNorm.get(i));
            popSizeParameter.setParameterValue(i,1.0);
//              popSizeParameter.setParameterValue(i,MathUtils.nextGaussian());
            }
    }

        public Parameter getPrecisionParameter() {
            return precisionParameter;
        }

        public Parameter getPopSizeParameter() {
            return popSizeParameter;
        }

    public Parameter getPopValue() {
        return popValue;
    }
    public Parameter getNumPoints() {
        return numPoints;
    }

    public Parameter getLambdaParameter() {
            return lambdaParameter;
        }

        public Parameter getLambdaBoundParameter() {
               return lambda_boundParameter;
        }

        public Parameter getChangePoints() {
            return changePoints;
        }

        public double getAlphaParameter(){
            return alphaParameter.getParameterValue(0);
        }

        public double getBetaParameter(){
            return betaParameter.getParameterValue(0);
        }


        public double [] getGPcoalfactor(){
            return GPcoalfactor;
        }


        public Parameter getcoalfactor(){
            return coalfactor;
        }

         public Parameter getCoalCounts(){
            return CoalCounts;
        }


        public Parameter getGPtype(){
            return GPtype;
        }


        public Parameter getGPcounts(){
            return GPcounts;
        }


	    public SymmTridiagMatrix getWeightMatrix() {
		    return weightMatrix.copy();
	}

//    Methods needed for GP-based

    public double [] getGPCoalInterval(){
        return GPCoalInterval;
    }

    public double [] getCoalTime(){
        return CoalTime;
    }

    public double getGPCoalInterval(int j){
        return GPCoalInterval[j];
    }

    public int [] getCoalPosIndicator() {
        return CoalPosIndicator;
    }

}

