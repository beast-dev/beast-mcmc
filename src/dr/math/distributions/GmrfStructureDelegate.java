/*
 * GaussianMarkovRandomField.java
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

package dr.math.distributions;

import dr.inference.distribution.RandomField;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 * @author Pratyusa Data
 */
public interface GmrfStructureDelegate  {

    GaussianMarkovRandomField.SymmetricTriDiagonalMatrix getQ();

//    void computeQ(GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q);

    abstract class Abstract extends AbstractModel implements GmrfStructureDelegate {

        final int dim;

        boolean qKnown;
        boolean savedQKnown;

        GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q;

        public Abstract(String name, int dim, Parameter lambda) {
            super(name);

            this.dim = dim;

            if (lambda != null) {
                addVariable(lambda);
            }
        }

        abstract void computeQ(GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q);

        @Override
        public GaussianMarkovRandomField.SymmetricTriDiagonalMatrix getQ() {

            if (!qKnown) {
                computeQ(q);
            }
            return q;
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            throw new RuntimeException("Not implemented");
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

    class Basic extends Abstract {

        public Basic(String name, int dim, Parameter lambda) {
            super(name, dim, lambda);
        }

        void computeQ(GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q) {

            q.diagonal[0] = 1.0;
            for (int i = 1; i < dim - 1; ++i) {
                q.diagonal[i] = 2.0;
            }
            q.diagonal[dim - 1] = 1.0;
            for (int i = 0; i < dim - 1; ++i) {
                q.offDiagonal[i] = -1.0;
            }
        }
    }

    class Weighted extends Basic {

        private final RandomField.WeightProvider weights;

        public Weighted(String name, int dim, Parameter lambda, RandomField.WeightProvider weights) {
            super(name, dim, lambda);
            this.weights = weights;
        }

        void computeQ(GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q) {
            q.diagonal[0] = weights.weight(0, 1);
            for (int i = 1; i < dim - 1; ++i) {
                q.diagonal[i] = weights.weight(i - 1, i) + weights.weight(i, i + 1);
            }
            q.diagonal[dim - 1] = weights.weight(dim - 2, dim - 1);
            for (int i = 0; i < dim - 1; ++i) {
                q.offDiagonal[i] = -weights.weight(i, i + 1);
            }
        }
    }

    class RhoConditionalAutoRegressive implements GmrfStructureDelegate {
        @Override
        public GaussianMarkovRandomField.SymmetricTriDiagonalMatrix getQ() {
            return null;
        }
    }

    class LSomethingConditionalAutoRegressive extends Abstract {

        private final GmrfStructureDelegate gmrf;

        public LSomethingConditionalAutoRegressive(String name,
                                                   GmrfStructureDelegate gmrf) {
            super(name, -1, null);
            this.gmrf = gmrf;
        }

        @Override
        void computeQ(GaussianMarkovRandomField.SymmetricTriDiagonalMatrix q) {

//            gmrf.computeQ(q);
//            GaussianMarkovRandomField.SymmetricTriDiagonalMatrix
//                    q = gmrf.getQ();
            // c * Q(1) + (1 - c) * I
//            return q;
        }
    }
}
