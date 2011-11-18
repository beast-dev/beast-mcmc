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

import com.sun.xml.internal.rngom.digested.DDataPattern;
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

import java.util.ArrayList;
import java.util.List;


/**
 * @author Vladimir Minin
 * @author Marc Suchard
 * @author Julia Palacios
 * @author Mandev
 */
public class GaussianProcessSkytrackLikelihood extends OldAbstractCoalescentLikelihood {

    protected Parameter popSizeParameter;
    protected Parameter groupSizeParameter;
    protected Parameter precisionParameter;
    protected Parameter lambda_boundParameter;
    protected Parameter lambdaParameter;
    protected Parameter betaParameter;

//    protected int GPfieldLength;
    protected double[] coalescentIntervals;
    protected double[] allintervals;
    protected double[] storedCoalescentIntervals;
    protected double[] sufficientStatistics;
    protected double[] allsufficient;
    protected double[] storedSufficientStatistics;

    protected double[] latentLocations;
    protected double[] storedlatentLocations;
    protected double[] latentpopSizeValues;
    protected double[] storedlatentpopSizeValues;

    protected double logFieldLikelihood;
    protected double storedLogFieldLikelihood;
    protected SymmTridiagMatrix weightMatrix;
	protected SymmTridiagMatrix storedWeightMatrix;
	protected MatrixParameter dMatrix;
	protected boolean rescaleByRootHeight;




    public GaussianProcessSkytrackLikelihood(Tree tree, Parameter popParameter, Parameter groupParameter,
                                             Parameter precParameter, Parameter lambda, Parameter beta,
                                             MatrixParameter dMatrix, /*boolean timeAwareSmoothing,*/
                                             boolean rescaleByRootHeight, /*Parameter latentPoints,*/ Parameter lambda_bound) {
        this(wrapTree(tree), popParameter, groupParameter, precParameter, lambda, beta, dMatrix, rescaleByRootHeight, lambda_bound);
    }


    public GaussianProcessSkytrackLikelihood(String name) {
		super(name);
	}

     private static List<Tree> wrapTree(Tree tree) {
        List<Tree> treeList = new ArrayList<Tree>();
        treeList.add(tree);
        return treeList;
    }

    public GaussianProcessSkytrackLikelihood(List<Tree> treeList, Parameter popParameter, Parameter groupParameter,
                                             Parameter precParameter, Parameter lambda, Parameter beta,
                                             MatrixParameter dMatrix,
                                             boolean rescaleByRootHeight, Parameter lambda_bound) {
        super(GaussianProcessSkytrackLikelihoodParser.SKYTRACK_LIKELIHOOD);



                this.popSizeParameter = popParameter;
                this.groupSizeParameter = groupParameter;
                this.precisionParameter = precParameter;
                this.lambdaParameter = lambda;
                this.betaParameter = beta;
                this.dMatrix = dMatrix;
                this.rescaleByRootHeight = rescaleByRootHeight;
                this.lambda_boundParameter= lambda_bound;

                addVariable(popSizeParameter);
                addVariable(precisionParameter);
                addVariable(lambdaParameter);
                addVariable(lambda_boundParameter);
                if (betaParameter != null) {
                    addVariable(betaParameter);
                }

                setTree(treeList);

                wrapSetupIntervals();

                int fieldLength = getCorrectFieldLength();
                int GPfieldLength= countChangePoints();

                allintervals = new double[GPfieldLength];
                allsufficient=new double[GPfieldLength];

                coalescentIntervals = new double[fieldLength];
		        storedCoalescentIntervals = new double[fieldLength];
		        sufficientStatistics = new double[fieldLength];
		        storedSufficientStatistics = new double[fieldLength];

                initializationReport();
                setupSufficientStatistics();


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
    }

    protected int getCorrectFieldLength() {
            return tree.getExternalNodeCount() - 1;
    }

    //
//
//         wrapSetupIntervals();
//         coalescentIntervals = new double[GPfieldLength];
//         storedCoalescentIntervals = new double[GPfieldLength];
//         sufficientStatistics = new double[GPfieldLength];
//         storedSufficientStatistics = new double[GPfieldLength];
//
////        setupGPWeights();
//
////        System.err.println("Aqui busco:"+GPfieldLength);
//

//
////    private void setupGPWeights() {
////        GPsetupSufficientStatistics();
////    }
//
//
//
//
//
//    public double calculateLogLikelihood() {
//        return 2.0;
// // TODO Return the correct log-density
//    }
//
//	public double getLogLikelihood() {
//		if (!likelihoodKnown) {
//			logLikelihood = calculateLogLikelihood();
//			likelihoodKnown = true;
//		}
//		return logLikelihood;
//	}
//
////    private final Parameter latentPoints;
//
//    private final Parameter lambda_bound;
//
//
//
    public void initializationReport() {

		System.out.println("Creating a GP based estimation of effective population trajectories:");
		System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
		System.out.println("\tIf you publish results using this model, please reference: Minin, Palacios, Suchard (XXXX), AAA");

	}
//
////    protected void setupGPWeights() {
////
////            setupSufficientStatistics();
////
////            //Set up the weight Matrix
////            double[] offdiag = new double[fieldLength - 1];
////            double[] diag = new double[fieldLength];
////
////            //First set up the offdiagonal entries;
////
////            if (!timeAwareSmoothing) {
////                for (int i = 0; i < fieldLength - 1; i++) {
////                    offdiag[i] = -1.0;
////                }
////
////
////            } else {
////                for (int i = 0; i < fieldLength - 1; i++) {
////                    offdiag[i] = -2.0 / (coalescentIntervals[i] + coalescentIntervals[i + 1]) * getFieldScalar();
////                }
////            }
////
////            //Then set up the diagonal entries;
////            for (int i = 1; i < fieldLength - 1; i++)
////                diag[i] = -(offdiag[i] + offdiag[i - 1]);
////
////            //Take care of the endpoints
////            diag[0] = -offdiag[0];
////            diag[fieldLength - 1] = -offdiag[fieldLength - 2];
////
////            weightMatrix = new SymmTridiagMatrix(diag, offdiag);
////        }
//
//
//
     //Sufficient Statistics for GP - coal+sampling
    // We do not consider getLineageCount == 1, we combine it with the next time
    // We do not need to make a difference between coalescent and sampling points
    //CoalescentIntervals here actually contain change points
    //sufficientStatistics here actually contain (k choose 2)

    protected void setupSufficientStatistics() {
	    int index = 0;
        int index2 = 0;

		double length = 0;
//		double weight = 0;
        System.err.println("Prueba dimension"+coalescentIntervals.length);
		for (int i = 0; i < getIntervalCount(); i++) {

//            if (getLineageCount(i)< 2) {
//                length+=getInterval(i);
//
//             } else {
			length += getInterval(i);
              System.err.println("Aqui va: getInterval:"+i+" "+getInterval(i)+ "with length "+length+
                    " with Type"+getIntervalType(i)+" lineage count "+getLineageCount(i));
            	allintervals[index] = length;
				allsufficient[index] =getLineageCount(i)*(getLineageCount(i)-1) / 2.0;
				index++;
                if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
				coalescentIntervals[index2] = length;
				sufficientStatistics[index2] = getLineageCount(i)*(getLineageCount(i)-1) / 2.0;
				index2++;
			        }
//                }
        }
        System.exit(-1);
    }

//    protected int countChangePoints(){
//        int countCPoints = 0;
//        for (int i =0; i<getIntervalCount();i++){
//            if (getLineageCount(i)>1){
//            countCPoints++;
//            }
//
//        }
//        return countCPoints;
//
//    }
  protected int countChangePoints(){
    return getIntervalCount();
    }
//
}
