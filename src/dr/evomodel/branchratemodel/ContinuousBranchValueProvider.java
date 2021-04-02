/*
 * CountableBranchCategoryProvider.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Xiang Ji
 */
public interface ContinuousBranchValueProvider {

    double getBranchValue(final Tree tree, final NodeRef node);

    class MidPoint implements ContinuousBranchValueProvider {

        @Override
        public double getBranchValue(Tree tree, NodeRef node) {

            assert (tree.getRoot() != node);

            final double midPoint = (tree.getNodeHeight(tree.getParent(node)) + tree.getNodeHeight(node)) * 0.5;

            return transform(midPoint);
        }

        double transform(double x) {
            return Math.log(x);
        }
    }

    class ConstrainedMidPoint extends AbstractModel implements ContinuousBranchValueProvider {

        final static String NAME = "constrainedMidPoint";
        private final Parameter heightLowerBound;

        public ConstrainedMidPoint(Parameter heightLowerBound) {
            super(NAME);
            this.heightLowerBound = heightLowerBound;
            addVariable(heightLowerBound);
        }

        @Override
        public double getBranchValue(Tree tree, NodeRef node) {
            assert (tree.getRoot() != node);

            final double midPoint;
            if (tree.getNodeHeight(tree.getParent(node)) > heightLowerBound.getParameterValue(0)) {
                midPoint = (tree.getNodeHeight(tree.getParent(node)) + tree.getNodeHeight(node)) * 0.5;
            } else {
                midPoint = heightLowerBound.getParameterValue(0);
            }

            return transform(midPoint);

        }

        double transform(double x) {
            return Math.log(x);
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {

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

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        }
    }
}
