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
import dr.evolution.tree.Tree;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import no.uib.cipr.matrix.SymmTridiagMatrix;

import java.util.List;

/**
 * @author Vladimir Minin
 * @author Marc Suchard
 * @author Julia Palacios
 * @author Mandev
 */
public class GaussianProcessSkytrackLikelihood extends GMRFSkyrideLikelihood {

    public GaussianProcessSkytrackLikelihood(List<Tree> treeList, Parameter popParameter, Parameter groupParameter,
                                             Parameter precParameter, Parameter lambda, Parameter beta,
                                             MatrixParameter dMatrix, /*boolean timeAwareSmoothing,*/
                                             boolean rescaleByRootHeight, /*Parameter latentPoints,*/ Parameter lambda_bound) {
        super(treeList, popParameter, groupParameter, precParameter, lambda, beta, dMatrix, /*timeAwareSmoothing,*/ false,
                rescaleByRootHeight);
//        this.latentPoints = latentPoints;
        this.lambda_bound = lambda_bound;
        addVariable(lambda_bound);
//        latentPoints.setDimension(1);
        initializationReport();


    }

//    protected double[] latentLocations;
//    protected double[] storedlatentLocations;
//    protected double[] latentpopSizeValues;
//    protected double[] storedlatentpopSizeValues;

    public double calculateLogLikelihood() {
        return 2.0; // TODO Return the correct log-density
    }

	public double getLogLikelihood() {
		if (!likelihoodKnown) {
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}
		return logLikelihood;
	}

//    private final Parameter latentPoints;

    private final Parameter lambda_bound;


	public void initializationReport() {
		System.out.println("Creating a GP based estimation of effective population trajectories:");
		System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
		System.out.println("\tIf you publish results using this model, please reference: Minin, Palacios, Suchard (XXXX), AAA");
	}

//    protected void setupGPWeights() {
//
//            setupSufficientStatistics();
//
//            //Set up the weight Matrix
//            double[] offdiag = new double[fieldLength - 1];
//            double[] diag = new double[fieldLength];
//
//            //First set up the offdiagonal entries;
//
//            if (!timeAwareSmoothing) {
//                for (int i = 0; i < fieldLength - 1; i++) {
//                    offdiag[i] = -1.0;
//                }
//
//
//            } else {
//                for (int i = 0; i < fieldLength - 1; i++) {
//                    offdiag[i] = -2.0 / (coalescentIntervals[i] + coalescentIntervals[i + 1]) * getFieldScalar();
//                }
//            }
//
//            //Then set up the diagonal entries;
//            for (int i = 1; i < fieldLength - 1; i++)
//                diag[i] = -(offdiag[i] + offdiag[i - 1]);
//
//            //Take care of the endpoints
//            diag[0] = -offdiag[0];
//            diag[fieldLength - 1] = -offdiag[fieldLength - 2];
//
//            weightMatrix = new SymmTridiagMatrix(diag, offdiag);
//        }


     //Sufficient Statistics for GP - coal+sampling
//    protected void setupSufficientStatistics() {
//	    int index = 0;
//
//		double length = 0;
//		double weight = 0;
//		for (int i = 0; i < getIntervalCount(); i++) {
//
//			length += getInterval(i);
//			weight += getInterval(i) * getLineageCount(i) * (getLineageCount(i) - 1);
//			if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
//				coalescentIntervals[index] = length;
//				sufficientStatistics[index] = weight / 2.0;
//				index++;
//				length = 0;
//				weight = 0;
//			}
//		}
//    }
//

}
