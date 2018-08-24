/*
 * ContinuousRateTransformation.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A. Suchard
 */

public interface ContinuousRateTransformation {

    double getNormalization();

    RateRescalingScheme getRateRescalingScheme();

    class Default extends AbstractModel implements ContinuousRateTransformation {
        private final Tree tree;

        private final RateRescalingScheme scheme;

        public Default(Tree tree, boolean scaleByTime, boolean useTreeLength) {
            super("ContinuousRateTransformation");
            this.tree = tree;

            scheme = scaleByTime ?
                    (useTreeLength ?
                            RateRescalingScheme.TREE_LENGTH :
                            RateRescalingScheme.TREE_HEIGHT) :
                    RateRescalingScheme.NONE;

        }

        private static final boolean DEBUG = false;

        public RateRescalingScheme getRateRescalingScheme() {
            return scheme;
        }

        @Override
        public double getNormalization() { // TODO Cache

            double norm = 1.0;
            switch (scheme) {
                case TREE_LENGTH:
                    norm = 1.0 / getTreeLength();
                    break;
                case TREE_HEIGHT:
                    norm = 1.0 / tree.getNodeHeight(tree.getRoot());
                    break;
            }

            if (DEBUG) {
                System.out.println("CRT.gN = " + norm);
            }

            return norm;
        }

        private double getTreeLength() {
            double treeLength = 0;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                if (!tree.isRoot(node)) {
                    treeLength += tree.getBranchLength(node);
                }
            }
            return treeLength;
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {

        }
        
        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        }

        @Override
        protected void storeState() {

        }

        @Override
        protected void restoreState() {

        }

        @Override
        protected void acceptState() {

        }
    }
}
