/*
 * HomogenousBranchSubstitutionModel.java
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
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BufferIndexHelper;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */

@Deprecated // Switching to BranchModel
public class HomogenousBranchSubstitutionModel extends AbstractModel implements BranchSubstitutionModel {
    public HomogenousBranchSubstitutionModel(SubstitutionModel substModel, FrequencyModel frequencyModel) {
        super("HomogenousBranchSubstitutionModel");

        this.substModel = substModel;
        addModel(substModel);
        this.frequencyModel = frequencyModel;
        addModel(frequencyModel);
    }

    /**
     * Homogenous model - returns the same substitution model for all branches/categories
     *
     * @param branchIndex
     * @param categoryIndex
     * @return
     */
    public EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex) {
        return substModel.getEigenDecomposition();
    }

    public void setEigenDecomposition(Beagle beagle, int eigenIndex, BufferIndexHelper bufferHelper, int dummy) {
        EigenDecomposition ed = getEigenDecomposition(eigenIndex, dummy);

        beagle.setEigenDecomposition(
                bufferHelper.getOffsetIndex(eigenIndex),
                ed.getEigenVectors(),
                ed.getInverseEigenVectors(),
                ed.getEigenValues());

    }

    public SubstitutionModel getSubstitutionModel(int branchIndex, int categoryIndex) {
        return substModel;
    }

    /**
     * Homogenous model - returns the same frequency model for all categories
     *
     * @param categoryIndex
     * @return
     */
    public double[] getStateFrequencies(int categoryIndex) {
        return frequencyModel.getFrequencies();
    }

    /**
     * Homogenous model - returns if substitution model can return complex diagonalization
     *
     * @return
     */
    public boolean canReturnComplexDiagonalization() {
        return substModel.canReturnComplexDiagonalization();
    }

    /**
     * Homogenous model - always returns model 0
     *
     * @param tree
     * @param node
     * @return
     */
    public int getBranchIndex(final Tree tree, final NodeRef node, int bufferIndex) {
        return 0;
    }

    public int getEigenCount() {
        return 1;
    }

    private final SubstitutionModel substModel;
    private final FrequencyModel frequencyModel;

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public void updateTransitionMatrices(Beagle beagle,
                                         int eigenIndex,
                                         BufferIndexHelper bufferHelper,
                                         final int[] probabilityIndices,
                                         final int[] firstDerivativeIndices,
                                         final int[] secondDervativeIndices,
                                         final double[] edgeLengths,
                                         int count) {
        beagle.updateTransitionMatrices(bufferHelper.getOffsetIndex(eigenIndex), probabilityIndices, firstDerivativeIndices,
                secondDervativeIndices, edgeLengths, count);

        //////////////////////////////////////////////////////

//		for(int k =0;k<probabilityIndices.length;k++){
//		
//		double tmp[] = new double[4 * 4 * 4];
//		beagle.getTransitionMatrix(probabilityIndices[k], // matrixIndex
//				tmp // outMatrix
//				);
//		
//		System.out.println(probabilityIndices[k]);
//		EpochBranchSubstitutionModel.printMatrix(tmp, 4, 4);
//		}

        //////////////////////////////////////////////////////

    }


    public int getExtraBufferCount(TreeModel treeModel) {
        // TODO Auto-generated method stub
        return 0;
    }


    public void setFirstBuffer(int bufferCount) {
        // TODO Auto-generated method stub

    }
}
