/*
 * EstimableStemWeightBranchSpecificBranchModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.branchmodel;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;
import java.util.logging.Logger;

/**
 * An extension of BranchSpecificBranchModel that allows for the stemWeight to be an estimated parameter.
 *
 * @author Andy Magee
 */
public class EstimableStemWeightBranchSpecificBranchModel extends BranchSpecificBranchModel {

    final private List<Parameter> stemWeightParameters = new ArrayList<>();
    protected Map<BitSet, Integer> stemWeightMap = new HashMap<>();
    private boolean hasBackbone = false;

    public EstimableStemWeightBranchSpecificBranchModel(TreeModel treeModel, SubstitutionModel rootSubstitutionModel) {
        super(treeModel, rootSubstitutionModel);
    }

    public void checkValidityOfBranchAssignments() {
        if (hasBackbone) {
            throw new RuntimeException("Not implemented yet");
        }

        if ( clades.size() > 1 ) {
            for (BitSet clade1 : clades.keySet()) {
                for (BitSet clade2 : clades.keySet()) {
                    if ( !clade1.equals(clade2) && clade1.intersects(clade2)) {
                        throw new RuntimeException("Overlapping clades are not allowed.");
                    }
                }
            }
        }

        Map<NodeRef, Mapping> externalNodeMap = getExternalNodeMap();
        if ( clades.size() > 0 && externalNodeMap.size() > 0 ) {
            for ( NodeRef node : externalNodeMap.keySet() ) {
                int nodeNumber = node.getNumber();
                for (BitSet clade : clades.keySet()) {
                    if ( clade.get(nodeNumber) ) {
                        throw new RuntimeException("Overlapping tip sets and clades are not allowed.");
                    }
                }
            }
        }

    }

    /**
     * Adds a substitution model specific to a clade.
     * @param taxonList a list of taxa who's MRCA define the clade
     * @param substitutionModel the substitution model
     * @param stemWeightParameter the proportion of the stem branch to include in this model (0, 1)
     * @throws TreeUtils.MissingTaxonException
     */
    public void addClade(TaxonList taxonList, SubstitutionModel substitutionModel, Parameter stemWeightParameter) throws TreeUtils.MissingTaxonException {
        List<SubstitutionModel> substitutionModels = getSubstitutionModels();
        int index = substitutionModels.indexOf(substitutionModel);
        if (index == -1) {
            index = substitutionModels.size();
            substitutionModels.add(substitutionModel);
            addModel(substitutionModel);
        }

        BitSet tips = TreeUtils.getTipsBitSetForTaxa(getTreeModel(), taxonList);
        Clade clade = new Clade(index, tips, stemWeightParameter.getParameterValue(0));
        clades.put(tips, clade);

        setRequiresMatrixConvolution(true);

        int stemIndex = stemWeightParameters.indexOf(stemWeightParameter);
        if (stemIndex == -1) {
            stemIndex = stemWeightParameters.size();
            stemWeightParameters.add(stemWeightParameter);
            addVariable(stemWeightParameter);
            stemWeightMap.put(tips, stemIndex);
        }

        Logger.getLogger("dr.evomodel.branchmodel")
                .info("\tAdding substitution model for clade defined by " + taxonList.getId() +
                        " with variable stem-weight parameter " + stemWeightParameter.getParameterName());
    }

    public void addBackbone(TaxonList taxonList, SubstitutionModel substitutionModel) throws TreeUtils.MissingTaxonException {
        hasBackbone = true;
        throw new UnsupportedOperationException("Not implemented yet");
    }


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable instanceof Parameter && stemWeightParameters.contains((Parameter)variable) && clades.size() > 0) {
            setUpdateNodeMaps(true);
        }
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void acceptState() {
    }

    void setupNodeMaps() {
        if (clades.size() > 0) {
            clearNodeMaps();
            updateStemWeights();
            super.setupNodeMaps();
        }
        setUpdateNodeMaps(false);
    }

    void processConvolvedBranch(NodeRef node, int index, double weight) {
        setConvolvedNodeMap(node, index, 0, weight);
    }

    void updateStemWeights() {
        double[] stemWeights = new double[stemWeightParameters.size()];
        for (int i = 0; i < stemWeightParameters.size(); i++) {
            stemWeights[i] = stemWeightParameters.get(i).getParameterValue(0);
        }
        for (BitSet clade : clades.keySet()) {
            clades.get(clade).setStemWeight(stemWeights[stemWeightMap.get(clade)]);
        }
    }
}