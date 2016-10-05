/*
 * DataSimulationDelegate.java
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

import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Variable;

import java.util.Collection;
import java.util.List;

/**
 * ProcessSimulationDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface ProcessSimulationDelegate extends ProcessOnTreeDelegate, Model, TreeTraitProvider {

    void simulate(List<DataLikelihoodDelegate.BranchOperation> branchOperations,
                  List<DataLikelihoodDelegate.NodeOperation> nodeOperations,
                  int rootNodeNumber);


    class Test implements ProcessSimulationDelegate {

        private final TreeTraversal treeTraversalDelegate;
        private TreeDataLikelihood treeDataLikelihood;
        private final Helper treeTraitProvider = new TreeTraitProvider.Helper();

        public Test(TreeModel treeModel, boolean[] updateNode, BranchRateModel branchRateModel) { // TODO deprecate
            this(new TreeTraversal(treeModel, branchRateModel, TreeTraversal.TraversalType.REVERSE_LEVEL_ORDER));
        }

        public Test(TreeTraversal treeTraversalDelegate) {
            this.treeTraversalDelegate = treeTraversalDelegate;
        }

        @Override
        public TreeTraversal.TraversalType getOptimalTraversalType() {
            return null;
        }

        @Override
        public void makeDirty() {

        }

        @Override
        public void storeState() {

        }

        @Override
        public void restoreState() {

        }

        @Override
        public void setCallback(TreeDataLikelihood treeDataLikelihood) {

        }

        /**
         * @return the id as a string.
         */
        @Override
        public String getId() {
            return null;
        }

        /**
         * set the id as a string.
         *
         * @param id
         */
        @Override
        public void setId(String id) {

        }

        /**
         * Adds a listener that is notified when the this model changes.
         *
         * @param listener
         */
        @Override
        public void addModelListener(ModelListener listener) {

        }

        /**
         * Remove a listener previously addeed by addModelListener
         *
         * @param listener
         */
        @Override
        public void removeModelListener(ModelListener listener) {

        }

        /**
         * This function should be called to store the state of the
         * entire model. This makes the model state invalid until either
         * an acceptModelState or restoreModelState is called.
         */
        @Override
        public void storeModelState() {

        }

        /**
         * This function should be called to restore the state of the entire model.
         */
        @Override
        public void restoreModelState() {

        }

        /**
         * This function should be called to accept the state of the entire model
         */
        @Override
        public void acceptModelState() {

        }

        /**
         * @return whether this model is in a valid state
         */
        @Override
        public boolean isValidState() {
            return false;
        }

        /**
         * @return the total number of sub-models
         */
        @Override
        public int getModelCount() {
            return 0;
        }

        /**
         * @param i
         * @return the ith sub-model
         */
        @Override
        public Model getModel(int i) {
            return null;
        }

        /**
         * @return the total number of variable in this model
         */
        @Override
        public int getVariableCount() {
            return 0;
        }

        /**
         * @param i
         * @return the ith variable
         */
        @Override
        public Variable getVariable(int i) {
            return null;
        }

        /**
         * @return the name of this model
         */
        @Override
        public String getModelName() {
            return null;
        }

        /**
         * is the model being listened to by another or by a likelihood?
         *
         * @return
         */
        @Override
        public boolean isUsed() {
            return false;
        }

        @Override
        public TreeTrait[] getTreeTraits() {
            return treeTraitProvider.getTreeTraits();
        }

        @Override
        public TreeTrait getTreeTrait(String key) {
            return treeTraitProvider.getTreeTrait(key);
        }

        @Override
        public void simulate(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations, int rootNodeNumber) {

        }
    }
}
