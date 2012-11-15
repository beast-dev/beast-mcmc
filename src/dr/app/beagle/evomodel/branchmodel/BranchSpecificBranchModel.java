/*
 * BranchSpecificBranchModel.java
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

package dr.app.beagle.evomodel.branchmodel;

import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.LocalClockModelParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * A branch model which allows different clades (defined by MRCAs of taxon lists) to have different
 * substitution models.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class BranchSpecificBranchModel extends AbstractModel implements BranchModel {

    private TreeModel treeModel;
    protected Map<BitSet, Clade> clades = new HashMap<BitSet, Clade>();

    private boolean updateNodeMaps = true;
    private Map<NodeRef, Mapping> nodeMap = new HashMap<NodeRef, Mapping>();

    private final SubstitutionModel rootSubstitutionModel;
    private final List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();

    private boolean requiresMatrixConvolution = false;

    public BranchSpecificBranchModel(TreeModel treeModel, SubstitutionModel rootSubstitutionModel) {

        super(LocalClockModelParser.LOCAL_CLOCK_MODEL);
        this.treeModel = treeModel;

        addModel(treeModel);

        this.rootSubstitutionModel = rootSubstitutionModel;
        addModel(rootSubstitutionModel);

        substitutionModels.add(rootSubstitutionModel);
    }

    /**
     * Adds a substitution model specific to a clade.
     * @param taxonList a list of taxa who's MRCA define the clade
     * @param substitutionModel the substitution model
     * @param stemWeight the proportion of the stem branch to include in this model (0, 1)
     * @throws Tree.MissingTaxonException
     */
    public void addClade(TaxonList taxonList, SubstitutionModel substitutionModel, double stemWeight) throws Tree.MissingTaxonException {
        int index = substitutionModels.indexOf(substitutionModel);
        if (index == -1) {
            index = substitutionModels.size();
            substitutionModels.add(substitutionModel);
            addModel(substitutionModel);
        }

        BitSet tips = Tree.Utils.getTipsForTaxa(treeModel, taxonList);
        Clade clade = new Clade(index, tips, stemWeight);
        clades.put(tips, clade);

        if (stemWeight > 0.0 || stemWeight < 1.0) {
            requiresMatrixConvolution = true;
        }
    }

    public void addExternalBranches(TaxonList taxonList, SubstitutionModel substitutionModel, double weight) throws Tree.MissingTaxonException {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    public void addBackbone(TaxonList taxonList, SubstitutionModel substitutionModel) throws Tree.MissingTaxonException {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public Mapping getBranchModelMapping(NodeRef branch) {
        if (updateNodeMaps) {
            setupNodeMaps();
        }
        Mapping mapping = nodeMap.get(branch);
        if (mapping != null) {
            return mapping;
        }
        return BranchModel.DEFAULT;
    }

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        return substitutionModels;
    }

    @Override
    public SubstitutionModel getRootSubstitutionModel() {
        return rootSubstitutionModel;
    }

    @Override
    public boolean requiresMatrixConvolution() {
        return requiresMatrixConvolution;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            updateNodeMaps = true;
        }
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
        updateNodeMaps = true;
    }

    protected void acceptState() {
    }


    private void setupNodeMaps() {
        setupNodeMaps(treeModel, treeModel.getRoot(), new BitSet());
    }

    private void setupNodeMaps(Tree tree, NodeRef node, BitSet tips) {
        Clade clade;

        if (tree.isExternal(node)) {
            tips.set(node.getNumber());
            clade = null;
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                BitSet childTips = new BitSet();
                setupNodeMaps(tree, child, childTips);

                tips.or(childTips);
            }
            clade = clades.get(tips);
        }

        if (clade != null) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                setNodeMap(tree, child, clade);
            }

            final double weight = clade.getStemWeight();

            if (weight > 0.0) {
                final int ancestralIndex;
                final int index = clade.getIndex();

                Mapping ancestoralMapping = nodeMap.get(node);

                if (ancestoralMapping != null) {
                    ancestralIndex = ancestoralMapping.getOrder()[0];
                } else {
                    ancestralIndex = 0;
                }

                nodeMap.put(node, new Mapping() {
                    @Override
                    public int[] getOrder() {
                        return new int[] { index , ancestralIndex };
                    }

                    @Override
                    public double[] getWeights() {
                        return new double[] { weight, 1.0 - weight };
                    }
                });

            }
        }
    }

    private void setNodeMap(Tree tree, NodeRef node, final Clade clade) {

        if (!tree.isExternal(node)) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                setNodeMap(tree, child, clade);
            }
        }

        nodeMap.put(node, new Mapping() {
            @Override
            public int[] getOrder() {
                return new int[] { clade.getIndex() };
            }

            @Override
            public double[] getWeights() {
                return new double[] { 1.0 };
            }
        });
    }

    private class Clade {

        Clade(int index, BitSet tips, double stemWeight) {
            this.index = index;
            this.tips = tips;
            this.stemWeight = stemWeight;
        }

        public int getIndex() {
            return index;
        }

        public BitSet getTips() {
            return tips;
        }

        public double getStemWeight() {
            return stemWeight;
        }

        private final int index;
        private final BitSet tips;
        private final double stemWeight;
    }

}