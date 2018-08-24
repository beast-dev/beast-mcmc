/*
 * GaussianProcessSkytrackLikelihood.java
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

package dr.evomodel.coalescent;

//import com.lowagie.text.Paragraph;

//import com.sun.servicetag.SystemEnvironment;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GaussianProcessSkytrackLikelihoodParser;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
//import dr.util.ComparableDouble;
//import dr.util.HeapSort;
import no.uib.cipr.matrix.*;
//import sun.font.TrueTypeFont;

import java.util.ArrayList;
import java.util.List;

//import dr.evolution.tree.TreeTrait;



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
//    protected Parameter numGridPoints;
    protected Parameter lambdaParameter;    //prior for lambda_bound, will be used in operators only
    protected Parameter betaParameter;
    protected Parameter alphaParameter;
    protected Parameter GPtype;
    protected Parameter GPcounts;
    protected Parameter coalfactor;
    protected Parameter popSizeParameter;     //before called GPvalues
    protected Parameter changePoints;
    protected Parameter Tmrca;
//    protected Parameter popValue;
    protected Parameter CoalCounts;
    protected Parameter numPoints;

//    protected double [] GPchangePoints;
//    protected double [] storedGPchangePoints;
    protected double [] GPcoalfactor;
    protected double [] storedGPcoalfactor;
    protected double [] GPCoalInterval;
    protected double [] storedGPCoalInterval;
    protected double [] backupIntervals;


    //    protected double [] storedcoalfactor;
//    protected int [] GPcounts;   //It changes values, no need to storage
//    protected int [] storedGPcounts;
    protected int [] CoalPosIndicator;
    protected int [] storedCoalPosIndicator;
    protected double [] CoalTime;
    protected double [] storedCoalTime;
    protected int numintervals;
    protected int numcoalpoints;
    protected double constlik;
    protected double storedconstlik;


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
                                             boolean rescaleByRootHeight,  Parameter lambda_bound, Parameter lambda_parameter, Parameter popParameter, Parameter alpha_parameter, Parameter beta_parameter, Parameter change_points, Parameter GPtype, Parameter GPcounts, Parameter coalfactor, Parameter CoalCounts, Parameter numPoints, Parameter Tmrca) {
        this(wrapTree(tree),  precParameter, rescaleByRootHeight, lambda_bound, lambda_parameter, popParameter, alpha_parameter, beta_parameter, change_points,GPtype,GPcounts,coalfactor,CoalCounts, numPoints, Tmrca);

    }


    public GaussianProcessSkytrackLikelihood(String name) {
		super(name);
	}




    public GaussianProcessSkytrackLikelihood(List<Tree> treeList,
                                             Parameter precParameter,
                                              boolean rescaleByRootHeight, Parameter lambda_bound, Parameter lambda_parameter, Parameter popParameter, Parameter alpha_parameter, Parameter beta_parameter, Parameter change_points, Parameter GPtype, Parameter GPcounts, Parameter coalfactor, Parameter CoalCounts, Parameter numPoints, Parameter Tmrca) {
        super(GaussianProcessSkytrackLikelihoodParser.SKYTRACK_LIKELIHOOD);



                this.popSizeParameter = popParameter;
                this.Tmrca=Tmrca;
//                this.popValue=popValues;
                this.changePoints=change_points;
                this.numPoints=numPoints;
//                this.groupSizeParameter = groupParameter;
                this.precisionParameter = precParameter;
                this.lambdaParameter = lambda_parameter;
                this.betaParameter = beta_parameter;
                this.alphaParameter=alpha_parameter;
//                this.dMatrix = dMatrix;
                this.rescaleByRootHeight = rescaleByRootHeight;
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

//          intervalCount = the size for constant vectors



//                int fieldLength = getCorrectFieldLength();
                numintervals= getIntervalCount();
                numcoalpoints=getCorrectFieldLength();

        GPcoalfactor = new double[numintervals];
        backupIntervals=new double[numintervals];
        GPCoalInterval=new double[numcoalpoints];
        storedGPCoalInterval=new double[numcoalpoints];
        CoalPosIndicator= new int[numcoalpoints];
        storedCoalPosIndicator=new int[numcoalpoints];
        CoalTime=new double[numcoalpoints];
        storedCoalTime=new double[numcoalpoints];


        storedGPcoalfactor = new double[numintervals];
        GPcounts.setDimension(numintervals);
        CoalCounts.setDimension(numcoalpoints);


//        storedGPcounts= new int[numintervals];
        GPtype.setDimension(numcoalpoints);
        numPoints.setParameterValue(0,numcoalpoints);
//        storedGPtype = new int[numcoalpoints];
        popSizeParameter.setDimension(numcoalpoints);


// NEED TO MOVE PopValue ---delete popValue also
//        int gridpoint= (int) numGridPoints.getParameterValue(0);

//        popValue.setDimension(gridpoint);

//        System.err.println("sets dimension");
        changePoints.setDimension(numcoalpoints);
        coalfactor.setDimension(numcoalpoints);
//        storedcoalfactor= new double[numcoalpoints];




                initializationReport();
                setupSufficientStatistics();
                setupGPvalues();
//                System.err.println("initial GP likelihood +priors"+getLogLikelihood());
//                System.err.println("like"+intervalsKnown);
//


//              System.err.println(getLogLikelihood());


         }


// Methods that override existent methods

    private boolean flagForJulia = false;

    /**
     * Demonstration of how to mark stuff as dirty
     * @param model
     * @param object
     * @param index
     */

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index); // Call super, since it may do something important
        if (model == tree) {
            // treeModel has changed; treeModel calls pushTreeChangedEvent that ultimately gets passed to here
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent tce = (TreeChangedEvent) object;
                // tce tells much about what type of event happened.  In general, one does not care.
//                System.err.println("Change in tree detected, flag true");
                flagForJulia = true; // flag set, so lazy work can occur elsewhere.
            } else {
                if (object instanceof Parameter){
//                    System.err.println("changing heights");
                    flagForJulia=true;
                }
                else {
                throw new IllegalArgumentException("Not sure what type of model changed event occurred: " + object.getClass().toString());
                }
            }
        }
    }

    public LogColumn[] getColumns() {
        // Add more LogColumn to the array if there are more things to log
        return new LogColumn[]{
                new VariableLengthColumn("changePoints", changePoints), new VariableLengthColumn("Gvalues",popSizeParameter)
        };
    }

    private class VariableLengthColumn extends LogColumn.Abstract {

        private final Parameter param;

        public VariableLengthColumn(String label, Parameter param) {
            super(label);
            this.param = param;
        }

        protected String getFormattedValue() {
            return convertToDelimited(param.getParameterValues());
        }

        // TODO The following functionality is generic and should be moved somewhere else and made static

        private static final String OPEN = "{";
        private static final String CLOSE = "}";
        private static final String DELIMIT = ",";

        private String convertToDelimited(double[] x) {
            StringBuilder sb = new StringBuilder(OPEN);
            final int dim = x.length;
            for (int i = 0; i < dim; ++i) {
                sb.append(Double.toString(x[i]));
                if (i < dim - 1) {
                    sb.append(DELIMIT);
                }
            }
            sb.append(CLOSE);
            return sb.toString();
        }
    }

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
        intervalsKnown = true;
    }



    //This is actually the Augmented loglikelihood for fixed genealogy. For sequence data directly
//    this becomes the coalescent point process prior on an augmented "tree"
    public double calculateLogLikelihood(Parameter Gfunction, Parameter latentCounts, Parameter eventType, Parameter upper_Bound, double [] Gfactor) {
        double upperBound = upper_Bound.getParameterValue(0);
        logGPLikelihood=-upperBound*getConstlik();

        for (int i=0; i<latentCounts.getSize(); i++){
            if (Gfactor[i]>0) {
                if (latentCounts.getParameterValue(i)<0){System.err.println("WARNING");}
          logGPLikelihood+=latentCounts.getParameterValue(i)*Math.log(upperBound*Gfactor[i]);
            }
        }

        double[] currentGfunction = Gfunction.getParameterValues();
        for (int i=0; i<Gfunction.getSize();i++){
        logGPLikelihood+= -Math.log(1+Math.exp(-eventType.getParameterValue(i)*currentGfunction[i]));
        }
        return logGPLikelihood;
    }


    public double getConstlik(){
        return constlik;
    }


//For fixed genealogy this contains the Augmented likelihood, the GP prior and prior on a the upper bound
	public double getLogLikelihood() {
		if (!likelihoodKnown) {
            if(flagForJulia) {
                System.err.println("recalculating intervals and counts");
                wrapSetupIntervals();
                recomputeValues();
                flagForJulia=false;
            }
           logLikelihood =
              calculateLogLikelihood(popSizeParameter,GPcounts,GPtype,lambda_boundParameter,GPcoalfactor)+calculateLogGP()
                        +getLogPriorLambda(lambdaParameter.getParameterValue(0),0.01,lambda_boundParameter.getParameterValue(0));
			likelihoodKnown = true;
		}
		return logLikelihood;
	}

    protected SymmTridiagMatrix getQmatrix(double precision, double[] x ) {
        SymmTridiagMatrix res;
        double trick=0.00000000001;
        double[] offdiag = new double[x.length - 1];
        double[] diag = new double[x.length];



        for (int i = 0; i < x.length - 1; i++) {
            offdiag[i] = precision*(-1.0 / (x[i+1]-x[i]));

            if (i< x.length-2){
                diag[i+1]= -offdiag[i]+precision*(1.0/(x[i+2]-x[i+1])+trick);
            }
        }
//              Diffuse prior correction - intrinsic
        //Take care of the endpoints
        diag[0] = -offdiag[0]+precision*trick;

        diag[x.length - 1] = -offdiag[x.length - 2]+precision*(trick);
        res = new SymmTridiagMatrix(diag, offdiag);
        return res;
    }


    //Calculates prior on g function
    protected double calculateLogGP() {
           SymmTridiagMatrix currentQ = getQmatrix(precisionParameter.getParameterValue(0), changePoints.getParameterValues());
           double currentLike;
           DenseVector diagonal1 = new DenseVector(popSizeParameter.getSize());
           DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());
           currentQ.mult(currentGamma, diagonal1);
           currentLike = -0.5 * logGeneralizedDeterminant(currentQ) - 0.5 * currentGamma.dot(diagonal1) - 0.5 * (popSizeParameter.getSize() - 1) * LOG_TWO_TIMES_PI;
//          System.err.println("the GP prior on f"+currentLike+" determinant"+logGeneralizedDeterminant(currentQ)+" size"+changePoints.getSize());
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
       	}


protected void restoreState() {
		super.restoreState();
  		System.arraycopy(storedGPcoalfactor, 0, GPcoalfactor, 0, storedGPcoalfactor.length);
        System.arraycopy(storedCoalTime,0,CoalTime,0,storedCoalTime.length);
        System.arraycopy(storedGPCoalInterval,0,GPCoalInterval,0,storedGPCoalInterval.length);
        System.arraycopy(storedCoalPosIndicator,0,CoalPosIndicator,0,storedCoalPosIndicator.length);
        constlik=storedconstlik;
    }

    protected void storeState() {
		super.storeState();
	 	System.arraycopy(GPcoalfactor, 0, storedGPcoalfactor, 0, GPcoalfactor.length);
        System.arraycopy(CoalTime,0,storedCoalTime,0,CoalTime.length);
		System.arraycopy(GPCoalInterval, 0, storedGPCoalInterval, 0, GPCoalInterval.length);
        System.arraycopy(CoalPosIndicator, 0, storedCoalPosIndicator,0,CoalPosIndicator.length);
        storedconstlik=constlik;
	}

       public String toString() {
        return getId() + "(" + Double.toString(getLogLikelihood()) + ")";
    }

    public void initializationReport() {

		System.out.println("Creating a GP based estimation of effective population trajectories:");
		System.out.println("\tIf you publish results using this model, please reference: Minin, Palacios, Suchard (XXXX), AAA");

	}

    public static void checkTree(TreeModel treeModel) {

            // todo Should only be run if there exists a zero-length interval - I don't actually understand why this is here

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

    protected void recomputeValues() {
        double length = 0.0;
        double prevLength=0.0;
        double count=0.0;
        int countcoal = 0;
        int pointer=0;
        int pointer2=0;
        constlik= 0;
        for (int i = 0; i < getIntervalCount(); i++) {
            length += getInterval(i);
            count=0.0;
            for (int j=pointer;j<changePoints.getSize();j++){

                if (changePoints.getParameterValue(j)<=length){

                     pointer++;
                     count++;

                }
            }
            GPcounts.setParameterValue(i,count);
            GPcoalfactor[i] =getLineageCount(i)*(getLineageCount(i)-1.0) / 2.0;
            constlik+=GPcoalfactor[i]*getInterval(i);

            if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                CoalPosIndicator[countcoal]=i;
                count=0;
                for (int s=pointer2;s<changePoints.getSize();s++){
                if (changePoints.getParameterValue(s)<=length){
                    pointer2++;
                    count++;
                } else
                    s=changePoints.getSize();
                }
                CoalCounts.setParameterValue(countcoal,count-1);
                CoalTime[countcoal]=length;
                GPCoalInterval[countcoal]=length-prevLength;
                coalfactor.setParameterValue(countcoal,getLineageCount(i)*(getLineageCount(i)-1)/2.0);
                countcoal++;
                prevLength=length;
            }
        }

//        TODO:This code should not be run in production, it is just for testing
        int sumcoal=0;
        int sumlat=0;

        for (int j=0;j<changePoints.getSize();j++){
            if( GPtype.getParameterValue(j)==1) {sumcoal++;}

        }
        for (int j=0;j<CoalCounts.getSize();j++){
            sumlat+=CoalCounts.getParameterValue(j);
        }
        if (sumcoal!=CoalCounts.getSize()){ System.err.println("WARNING CONSISTENCY 1"); }
        if (sumlat!=(changePoints.getSize()-CoalCounts.getSize())){ System.err.println("WARNING CONSISTENCY 2:"+sumlat+"and changePts size"+changePoints.getSize());}
        Tmrca.setParameterValue(0,CoalTime[countcoal-1]);
    }

    protected void setupSufficientStatistics() {
		double length = 0.0;
        double prevLength=0.0;
        int countcoal = 0;
        constlik= 0;
		for (int i = 0; i < getIntervalCount(); i++) {

			  length += getInterval(i);
              GPcounts.setParameterValue(i,0.0);
              GPcoalfactor[i] =getLineageCount(i)*(getLineageCount(i)-1.0) / 2.0;
              constlik+=GPcoalfactor[i]*getInterval(i);
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
        Tmrca.setParameterValue(0,CoalTime[countcoal-1]);

    }


    protected int getCorrectFieldLength() {
        return tree.getExternalNodeCount() - 1;
    }

    protected void setupQmatrix(double precision) {
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
        }
        UpperSPDBandMatrix Qcurrent = new UpperSPDBandMatrix(weightMatrix, 1);
        BandCholesky U = new BandCholesky(length,1,true);
        U.factor(Qcurrent);
        UpperTriangBandMatrix CholeskyUpper = U.getU();
        CholeskyUpper.solve(StandNorm,MultiNorm);
        for (int i=0; i<length;i++){
            popSizeParameter.setParameterValue(i,1.0);
        }
    }

        public Parameter getPrecisionParameter() {
            return precisionParameter;
        }

        public Parameter getPopSizeParameter() {
            return popSizeParameter;
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

