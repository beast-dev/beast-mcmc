/*
 * CountableBranchCategoryProvider.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeDoubleTraitProvider;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public interface CountableBranchCategoryProvider extends TreeTrait<Double> {

    public int getBranchCategory(final Tree tree, final NodeRef node);

    public void setCategoryCount(final int count);

    public int getCategoryCount();

    public class SingleBranchCategoryModel implements CountableBranchCategoryProvider {
        @Override
        public int getBranchCategory(final Tree tree, final NodeRef node) {
            return 0;
        }
        @Override
        public void setCategoryCount(final int count) {
            // Do nothing
        }
        @Override
        public int getCategoryCount() {
            return 1;
        }

        @Override
        public String getTraitName() {
            return "categories";
        }

        @Override
        public Intent getIntent() {
            return Intent.BRANCH;
        }

        @Override
        public Class getTraitClass() {
            return Integer.class;
        }

        @Override
        public Double getTrait(Tree tree, NodeRef node) {
            return 1.0;
        }

        @Override
        public String getTraitString(Tree tree, NodeRef node) {
            return "1";
        }

        @Override
        public boolean getLoggable() {
            return true;
        }
    }

    public abstract class BranchCategoryModel extends TreeParameterModel implements CountableBranchCategoryProvider {

        public BranchCategoryModel(TreeModel tree, Parameter parameter) {
            super(tree, parameter, false, Intent.BRANCH);

            this.categoryParameter = parameter;
//            for (int i = 0; i < parameter.getDimension(); ++i) {
//                categoryParameter.setParameterValue(i, 0.0);
//            }
            this.treeModel = tree;

            categoryCount = 1;
        }

//        public BranchCategoryModel(TreeModel tree, Parameter parameter, boolean resetCategories) {
//            super(tree, parameter, false, Intent.BRANCH);
//
//            this.categoryParameter = parameter;
//            this.treeModel = tree;
//        }

		public void setCategoryCount(final int count) {

            categoryCount = count;
            Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, categoryParameter.getDimension());
            categoryParameter.addBounds(bound);

            for (int i = 0; i < categoryParameter.getDimension(); ++i) {
                if (categoryParameter.getParameterValue(i) >= categoryCount) {
                    categoryParameter.setParameterValue(i, categoryCount - 1);
                }
            }
        }

        @Override
        public int getBranchCategory(final Tree tree, final NodeRef node) {
            return (int) Math.round(getNodeValue(tree, node));
        }

        @Override
        public int getCategoryCount() {
            return categoryCount;
        }

        protected int categoryCount;
        protected final Parameter categoryParameter;
        protected final TreeModel treeModel;
    }

    public class IndependentBranchCategoryModel extends BranchCategoryModel {

        public IndependentBranchCategoryModel(TreeModel tree, Parameter parameter) {
            super(tree, parameter);
        }

        public void randomize() {
            for (NodeRef node : treeModel.getNodes()) {
                if (node != treeModel.getRoot()) {
                    int index = MathUtils.nextInt(categoryCount);
                    setNodeValue(treeModel, node, index);
                }
            }
        }
    }

    public class MarkovJumpBranchCategoryModel extends BranchCategoryModel {

        public MarkovJumpBranchCategoryModel(MarkovJumpsTraitProvider markovJumpTrait, Parameter parameter) {
            super(markovJumpTrait.getTreeModel(), parameter);
        }

        @Override
        public int getBranchCategory(final Tree tree, final NodeRef node) {
            synchronized (this) {
                if (traitsChanged) {
                    updateTraitRateCategories();
                    traitsChanged = false;
                }
            }
            return super.getBranchCategory(tree, node);
        }

        private void updateTraitRateCategories() {

        }

        public void handleModelChangedEvent(Model model, Object object, int index) {
            if (model == treeModel) {
                traitsChanged = true;
                fireModelChanged();
            } else {
                throw new IllegalArgumentException("Unknown model component!");
            }
        }

        private boolean traitsChanged = true;
    }

    public class CladeBranchCategoryModel extends BranchCategoryModel {

        public CladeBranchCategoryModel(TreeModel tree, Parameter parameter) {
            super(tree, parameter);
        }

//        public CladeBranchCategoryModel(TreeModel treeModel,
//				Parameter categories, boolean resetCategories) {
//        	super(treeModel, categories, resetCategories);
//		}

		public void handleModelChangedEvent(Model model, Object object, int index) {
            if (model == treeModel) {
                cladesChanged = true;
                fireModelChanged();
            } else {
                throw new IllegalArgumentException("Unknown model component!");
            }
        }

        private void updateCladeRateCategories() {
            if (leafSetList != null) {
                for (NodeRef node : treeModel.getNodes()) {
                    if (node != treeModel.getRoot()) {
                        setNodeValue(treeModel, node, 0.0);
                    }
                }
                for (CladeContainer clade : leafSetList) {
                    NodeRef node = Tree.Utils.getCommonAncestorNode(treeModel, clade.getLeafSet());
                    if (node != treeModel.getRoot()) {
                        setNodeValue(treeModel, node, clade.getRateCategory());
                    }
                }
            }

            if (trunkSetList != null) {
                //we keep the default rates assignments by clade definitions if they exist (leafSetList != null), if they do not exist, set default to 0.0
                if (leafSetList == null) {
                    for (NodeRef node : treeModel.getNodes()) {
                        if (node != treeModel.getRoot()) {
                            setNodeValue(treeModel, node, 0.0);
                        }
                    }
                }
                // currently, specific backbone rates will overwrite branch assignments by clade definitions
                //TODO: think about turning this around. One can imagine setting backbone rates and then additional rates based on clade definitions
                for (CladeContainer trunk : trunkSetList) {
                    for (NodeRef node : treeModel.getNodes()) {
                        if (onAncestralPath(treeModel, node, trunk.getLeafSet())) {
                            if (node != treeModel.getRoot()) {
                                setNodeValue(treeModel, node, trunk.getRateCategory());
                            }
                        }
                    }
                }
            }
        }

        private boolean onAncestralPath(Tree tree, NodeRef node, Set targetSet) {

            if (tree.isExternal(node)) return false;

            Set leafSet = Tree.Utils.getDescendantLeaves(tree, node);
            int size = leafSet.size();

            leafSet.retainAll(targetSet);

            if (leafSet.size() > 0) {

                // if all leaves below are in target then check just above.
                if (leafSet.size() == size) {

                    Set superLeafSet = Tree.Utils.getDescendantLeaves(tree, tree.getParent(node));
                    superLeafSet.removeAll(targetSet);

                    // the branch is on ancestral path if the super tree has some non-targets in it
                    return (superLeafSet.size() > 0);

                } else return true;

            } else return false;
        }

        public void setClade(TaxonList taxonList, int rateCategory, boolean includeStem, boolean excludeClade, boolean trunk) throws Tree.MissingTaxonException {
            if (!excludeClade) {
                throw new IllegalArgumentException("Including clades not yet implemented in countable branch rate mixture models.");
            }

            if (!includeStem) {
                throw new IllegalArgumentException("Excluding stems not yet implemented in countable branch rate mixture models.");
            }

            Set<String> leafSet = Tree.Utils.getLeavesForTaxa(treeModel, taxonList);
            if (!trunk) {
                if (leafSetList == null) {
                    leafSetList = new ArrayList<CladeContainer>();
                }
                leafSetList.add(new CladeContainer(leafSet, rateCategory));
                cladesChanged = true;
            } else {
                if (trunkSetList == null) {
                    trunkSetList = new ArrayList<CladeContainer>();
                }
                trunkSetList.add(new CladeContainer(leafSet, rateCategory));
                cladesChanged = true;
            }
            if (rateCategory >= categoryCount) {
                categoryCount = rateCategory + 1;
            }
        }


        @Override
        public int getBranchCategory(final Tree tree, final NodeRef node) {
            synchronized (this) {
                if (cladesChanged) {
                    updateCladeRateCategories();
                    cladesChanged = false;
                }
            }
            return super.getBranchCategory(tree, node);
        }

        private class CladeContainer {
            private Set<String> leafSet;
            private int rateCategory;

            public CladeContainer(Set<String> leafSet, int rateCategory) {
                this.leafSet = leafSet;
                this.rateCategory = rateCategory;
            }

            public Set<String> getLeafSet() {
                return leafSet;
            }

            public int getRateCategory() {
                return rateCategory;
            }
        }

        private boolean cladesChanged = false;
        private List<CladeContainer> leafSetList = null;
        private List<CladeContainer> trunkSetList = null;
    }
}
