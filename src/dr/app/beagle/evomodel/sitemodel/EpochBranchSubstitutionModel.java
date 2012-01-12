/*
 * EpochBranchSubstitutionModel.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Filip Bielejec
 * @author Marc A. Suchard
 * @version $Id$
 */
public class EpochBranchSubstitutionModel extends AbstractModel implements BranchSubstitutionModel, Citable {

    // TODO Most of the functions below are a copy-n-paste from ExternalInternalBranchSubstitutionModel and need
    // TODO to be rewritten.

    public EpochBranchSubstitutionModel(List<SubstitutionModel> substModelList, List<FrequencyModel> frequencyModelList,
                                        Parameter epochTimes) {
        super("EpochBranchSubstitutionModel");

        if (substModelList.size() != 2) {
            throw new IllegalArgumentException("EpochBranchSubstitutionModel requires two SubstitutionModels");
        }

        if (frequencyModelList.size() != 1) {
            throw new IllegalArgumentException("EpochBranchSubstitutionModel requires one FrequencyModel");
        }

        this.substModelList = substModelList;
        this.frequencyModelList = frequencyModelList;
        this.epochTimes = epochTimes;

        for (SubstitutionModel model : substModelList) {
            addModel(model);
        }
        for (FrequencyModel model : frequencyModelList) {
            addModel(model);
        }

        addVariable(epochTimes);
    }

    public int getBranchIndex(final Tree tree, final NodeRef node) {
        // TODO return substModelList.size() if branch will require convolution
        return 0;
    }

    public EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex) {
        return substModelList.get(branchIndex).getEigenDecomposition();
    }

    public SubstitutionModel getSubstitutionModel(int branchIndex, int categoryIndex) {
        return substModelList.get(branchIndex);
    }

    public double[] getStateFrequencies(int categoryIndex) {
        return frequencyModelList.get(categoryIndex).getFrequencies();
    }

    public boolean canReturnComplexDiagonalization() {
        for (SubstitutionModel model : substModelList) {
            if (model.canReturnComplexDiagonalization()) {
                return true;
            }
        }
        return false;
    }

    public int getEigenCount() {
        return substModelList.size() + 1; // Use an extra eigenIndex to identify branches that need convolution
    }

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

    public void updateTransitionMatrices(Beagle beagle, int eigenIndex, final int[] probabilityIndices,
                                         final int[] firstDerivativeIndices, final int[] secondDervativeIndices,
                                         final double[] edgeLengths, int count) {

        if (eigenIndex < substModelList.size()) {
            // Branch falls in a single category
            beagle.updateTransitionMatrices(eigenIndex, probabilityIndices, firstDerivativeIndices,
                    secondDervativeIndices, edgeLengths, count);
        } else {
            // Branch requires convolution of two or more matrices
            // TODO Do fun work here
        }
    }

    private final List<SubstitutionModel> substModelList;
    private final List<FrequencyModel> frequencyModelList;
    private final Parameter epochTimes;

    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("F", "Bielejec"),
                                new Author("P", "Lemey"),
                                new Author("MA", "Suchard")
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
        return citations;
    }
}