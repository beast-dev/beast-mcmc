/*
 * DataLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood;

import dr.inference.model.Model;
import dr.xml.Reportable;

import java.util.List;

/**
 * DataLikelihoodDelegate - interface for a plugin delegate for the data likelihood.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface DataLikelihoodDelegate extends ProcessOnTreeDelegate, Model, Reportable {

    void makeDirty();

    void storeState();

    void restoreState();

    double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations,
                               int rootNodeNumber) throws LikelihoodException;

    int getTraitCount();

    int getTraitDim();

    RateRescalingScheme getRateRescalingScheme();

    class LikelihoodException extends Exception { }

    class LikelihoodUnderflowException extends LikelihoodException { }

    class LikelihoodRescalingException extends LikelihoodException { }

    void setCallback(TreeDataLikelihood treeDataLikelihood);

    int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations);

//    int getActiveNodeIndex(final int index);
//
//    int getActiveMatrixIndex(final int index);
}
