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

//import com.sun.xml.internal.rngom.digested.DDataPattern;
import dr.app.beast.BeastDialog;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GaussianProcessSkytrackLikelihoodParser;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SymmTridiagEVD;
import no.uib.cipr.matrix.SymmTridiagMatrix;
import sun.jvm.hotspot.memory.SystemDictionary;


import java.util.ArrayList;
import java.util.List;
//import java.util.Hashtable;


/**
 * @author Vladimir Minin
 * @author Marc Suchard
 * @author Julia Palacios
 * @author Mandev
 */
public class GaussianProcessSkytrackLikelihood extends OldAbstractCoalescentLikelihood {

//    protected Parameter popSizeParameter;
//    protected Parameter groupSizeParameter;
    protected Parameter precisionParameter;
    protected Parameter lambda_boundParameter;
    protected Parameter numGridPoints;
    protected Parameter lambdaParameter;    //prior for lambda_bound, will be used in operators only
//    protected Parameter lambdaParameter;
//    protected Parameter betaParameter;

//   Those that do not change in size  - fixed per tree -hence need to store/restore
    protected double[] GPchangePoints; //s
    protected double [] storedGPchangePoints;
    protected double [] GPcoalfactor;
    protected double [] storedGPcoalfactor;
    protected int [] GPcounts;   //It changes values, no need to storage
//    protected int [] storedGPcounts;
    protected int numintervals;
    protected double constlik;

//    Those that change size, they are initialized per tree, no need to store them
//    use as Parameter since they will be changing by operators
//    protected Parameter GPtimepoints;  //tree + latent
//    protected double GPintervalkey;         // membership that links with those that do not change in size
//    protected Parameter GPcoalfactor2;        // choose(k,2) depending on membership
    protected int[] GPtype;  // 1 if observed, -1 if latent
    public double[] GPvalues;     //may need to change type: Parameter? didn't know how to work with it

    protected double logGPLikelihood;
    protected double storedLogGPLikelihood;
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
                                             boolean rescaleByRootHeight, Parameter numGridPoints,  Parameter lambda_bound) {
        this(wrapTree(tree),  precParameter, rescaleByRootHeight, numGridPoints, lambda_bound);
    }


    public GaussianProcessSkytrackLikelihood(String name) {
		super(name);
	}




    public GaussianProcessSkytrackLikelihood(List<Tree> treeList,
                                             Parameter precParameter,
                                              boolean rescaleByRootHeight, Parameter numGridPoints, Parameter lambda_bound) {
        super(GaussianProcessSkytrackLikelihoodParser.SKYTRACK_LIKELIHOOD);



//                this.popSizeParameter = popParameter;
//                this.groupSizeParameter = groupParameter;
                this.precisionParameter = precParameter;
//                this.lambdaParameter = lambda;
//                this.betaParameter = beta;
//                this.dMatrix = dMatrix;
                this.rescaleByRootHeight = rescaleByRootHeight;
                this.numGridPoints = numGridPoints;
                this.lambda_boundParameter= lambda_bound;

//                addVariable(GPvalues);
                addVariable(precisionParameter);
//                addVariable(lambdaParameter);
//                addVariable(lambda_boundParameter);
//                if (betaParameter != null) {
//                    addVariable(betaParameter);
//                }

                setTree(treeList);

                wrapSetupIntervals();
//        intervalCount = the size for constant vectors



//                int fieldLength = getCorrectFieldLength();
                int numintervals= getIntervalCount();

        GPchangePoints = new double[numintervals];
        storedGPchangePoints = new double[numintervals];
        GPcoalfactor = new double[numintervals];
        storedGPcoalfactor = new double[numintervals];
        GPcounts = new int[numintervals];
        GPtype=new int[numintervals];




                initializationReport();
                setupSufficientStatistics();
                setupGPvalues();  //need to be changed

//
//                weightMatrix = new SymmTridiagMatrix(numintervals);


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
// Is it ok to have this public and override the one in OldAbstract...
    //I will use specific input
    public double calculateLogLikelihood(double[] Gfunction, SymmTridiagMatrix QMatrix, int[] latentCounts, int [] eventType, double upperBound, double [] Gintervals, int [] Gfactor) {
        logGPLikelihood=-getConstlik();

        for (int i=0; i<Gfunction.length;i++){
        logGPLikelihood+= latentCounts[i]*Math.log(upperBound*Gfactor[i])-upperBound*Gfactor[i]*Gintervals[i]-Math.log(1+Math.exp(-eventType[i]*Gfunction[i])   );
        }
        return 0.0;
// TODO Return the correct log-density  (augmented?)
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
	public double getLogLikelihood() {
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}
		return logLikelihood;
//        return 0.0;
	}

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type){
		likelihoodKnown = false;
        // Parameters (precision and popsizes do not change intervals or GMRF Q matrix (I DON'T UNDERSTAND THIS)
	}


protected void restoreState() {
		super.restoreState();
		System.arraycopy(storedGPchangePoints, 0, GPchangePoints, 0, storedGPchangePoints.length);
		System.arraycopy(storedGPcoalfactor, 0, GPcoalfactor, 0, storedGPcoalfactor.length);
//		weightMatrix = storedWeightMatrix;
        logGPLikelihood = storedLogGPLikelihood;
    }



protected void storeState() {
		super.storeState();
		System.arraycopy(GPchangePoints, 0, storedGPchangePoints, 0, GPchangePoints.length);
		System.arraycopy(GPcoalfactor, 0, storedGPcoalfactor, 0, GPcoalfactor.length);
//		storedWeightMatrix = weightMatrix.copy();
        storedLogGPLikelihood = logGPLikelihood;
	}
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


		double length = 0;
        constlik= 0;
		for (int i = 0; i < getIntervalCount(); i++) {

			length += getInterval(i);
//              System.err.println(i+"  "+length+"  "+getLineageCount(i)+" type "+getIntervalType(i));
                GPcounts[i]=0;
                GPtype[i]=1;
                GPchangePoints[i]=length;
             	GPcoalfactor[i] =getLineageCount(i)*(getLineageCount(i)-1) / 2.0;
                constlik+=GPcoalfactor[i]*getInterval(i);
                if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                    GPcounts[i]=1;
			        }
        }


    }

    protected void setupQmatrix(double precision) {


            //Set up the weight Matrix
            double[] offdiag = new double[getIntervalCount() - 1];
            double[] diag = new double[getIntervalCount()];


             for (int i = 0; i < getIntervalCount() - 1; i++) {
                    offdiag[i] = precision*(-1.0 / getInterval(i));
                    diag[i+1]= precision*(1.0/getInterval(i)+1.0/getInterval(i+1)+.000001);

                }
//              Diffuse prior correction - intrinsic
             //Take care of the endpoints
            diag[0] = -offdiag[0]+precision*.000001;
            diag[getIntervalCount() - 1] = -offdiag[getIntervalCount() - 2]+precision*(.000001);

            weightMatrix = new SymmTridiagMatrix(diag, offdiag);
        }



    protected void setupGPvalues() {
        setupQmatrix(precisionParameter.getParameterValue(0));
        GPvalues = new double[numintervals];

        // Need to add L<-chol.spam(Q), y2<-rnorm(nrow(Q),0,1), g<-backsolve(L,y2)

    }




//    Methods needed for GP-based

}
