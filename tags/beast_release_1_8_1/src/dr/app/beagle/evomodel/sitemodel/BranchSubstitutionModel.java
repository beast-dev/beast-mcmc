/*
 * BranchSubstitutionModel.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.sitemodel;

import beagle.Beagle;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc A. Suchard
 * @version $Id$
 */

@Deprecated // Switching to BranchModel
public interface BranchSubstitutionModel extends Model {

    EigenDecomposition getEigenDecomposition(int modelIndex, int categoryIndex);

    SubstitutionModel getSubstitutionModel(int modelIndex, int categoryIndex);

    double[] getStateFrequencies(int categoryIndex);

    public int getBranchIndex(final Tree tree, final NodeRef node, int bufferIndex);

    public int getEigenCount();

    boolean canReturnComplexDiagonalization();

    void updateTransitionMatrices(
            Beagle beagle,
            int eigenIndex,
            BufferIndexHelper bufferHelper,
            final int[] probabilityIndices,
            final int[] firstDerivativeIndices,
            final int[] secondDervativeIndices,
            final double[] edgeLengths,
            int count);

	int getExtraBufferCount(TreeModel treeModel);

	void setFirstBuffer(int bufferCount);

	void setEigenDecomposition(Beagle beagle, int eigenIndex, BufferIndexHelper bufferHelper, int dummy);
	
	

}
