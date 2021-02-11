/*
 * RateProvider.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.hawkes;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Andrew Holbrook
 * @author Xiang Ji
 * @author Marc Suchard
 */
public interface HawkesRateProvider {

    void setRandomRates(HawkesCore hawkesCore);

    Parameter getParameter();

    double[] orderByNodeIndex(double[] orderByTime);

    class None implements HawkesRateProvider {

        @Override
        public void setRandomRates(HawkesCore hawkesCore) {
            // do nothing
        }

        @Override
        public Parameter getParameter() {
            throw new RuntimeException("No rate parameter in the 'None' case of Hawkes");
        }

        @Override
        public double[] orderByNodeIndex(double[] orderByTime) {
            throw new RuntimeException("No rate parameter in the 'None' case of Hawkes");
        }
    }

    class Default extends AbstractModel implements HawkesRateProvider {

        private Parameter rate;
        private int[] indices;
        private boolean[] onTree;

        Default(Parameter rate,
                int[] indices,
                boolean[] onTree) {
            super("HawkesRateProvider$Default");
            this.rate = rate;
            this.indices = indices;
            this.onTree = onTree;
            addVariable(rate);
        }

        @Override
        public void setRandomRates(HawkesCore hawkesCore) {
            double[] orderedRates = new double[rate.getDimension()];
            for (int i = 0; i < rate.getDimension(); i++) {
                orderedRates[i] = rate.getParameterValue(indices[i]);
            }
            hawkesCore.setRandomRates(orderedRates);
        }

        public double[] orderByNodeIndex(double[] orderByTime) {
            double[] nodeOrdered = new double[orderByTime.length];
            for (int i = 0; i < orderByTime.length; i++) {
                nodeOrdered[indices[i]] = orderByTime[i];
            }
            return nodeOrdered;
        }

        @Override
        public Parameter getParameter() {
            return rate;
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
