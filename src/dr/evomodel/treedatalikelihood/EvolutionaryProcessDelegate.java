/*
 * EvolutionaryProcessDelegate.java
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

import beagle.Beagle;
import dr.evomodel.substmodel.SubstitutionModel;

/**
 * Implementations of this interface are used to delegate control of substitution
 * models along branches of a tree.
 *
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface EvolutionaryProcessDelegate {

    boolean canReturnComplexDiagonalization();

    int getEigenBufferCount();

    int getMatrixBufferCount();

    int getSubstitutionModelCount();

    SubstitutionModel getSubstitutionModel(int index);

    int getEigenIndex(int bufferIndex);

    int getMatrixIndex(int branchIndex);

    double[] getRootStateFrequencies();

    void updateSubstitutionModels(Beagle beagle, boolean flipBuffers);

    void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLengths, int updateCount, boolean flipBuffers);

    void flipTransitionMatrices(int[] branchIndices, int updateCount);

    void storeState();

    void restoreState();
}
