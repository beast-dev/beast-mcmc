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

import dr.evomodel.branchratemodel.ContinuousBranchValueProvider;
import dr.evomodel.tree.TreeModel;
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

    void updateRateGradient(double[] gradient);

    void updateRateHessian(double[] gradient, double[] hessian, double[] result);

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

        @Override
        public void updateRateGradient(double[] gradient) {
            // do nothing
        }

        @Override
        public void updateRateHessian(double[] gradient, double[] hessian, double[] result) {
            // do nothing
        }
    }

    class Default extends AbstractModel implements HawkesRateProvider {

        protected Parameter rate;
        protected int[] indices;
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
            hawkesCore.setRandomRates(orderByTime(getRatesByNodes()));
        }

        protected double[] getRatesByNodes() {
            return rate.getParameterValues();
        }

        public double[] orderByTime(double[] orderByNodes) {
            double[] timeOrdered = new double[orderByNodes.length];
            for (int i = 0; i < orderByNodes.length; i++) {
                timeOrdered[i] = orderByNodes[indices[i]];
            }
            return timeOrdered;
        }

        public double[] orderByNodeIndex(double[] orderByTime) {
            double[] nodeOrdered = new double[orderByTime.length];
            for (int i = 0; i < orderByTime.length; i++) {
                nodeOrdered[indices[i]] = orderByTime[i];
            }
            return nodeOrdered;
        }

        @Override
        public void updateRateGradient(double[] gradient) {
            // do nothing
        }

        @Override
        public void updateRateHessian(double[] gradient, double[] hessian, double[] result) {
            // do nothing
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

    class GLM extends Default {

        private Parameter coefficients;
        private TreeModel tree;
        private double[] nodeTimes;
        private boolean timeEffect;
        private boolean intercept;
        private ContinuousBranchValueProvider branchValueProvider;
        private final Options options;

        GLM(Parameter rate, Parameter coefficients, TreeModel tree, int[] indices, double[] orderedTimes,
            boolean[] onTree, boolean timeEffect, boolean intercept) {

            super(rate, indices, onTree);

            this.coefficients = coefficients;
            this.tree = tree;
            this.timeEffect = timeEffect;
            this.intercept = intercept;
            this.nodeTimes = orderByNodeIndex(orderedTimes);
            this.branchValueProvider = new ContinuousBranchValueProvider.MidPoint();
            this.options = new Options(intercept, timeEffect);

            addVariable(coefficients);
        }

        protected double[] getRatesByNodes() {
            double[] residues = super.getRatesByNodes();
            double[] rates = new double[rate.getDimension()];
            for (int i = 0; i < rate.getDimension(); i++) {
                rates[i] = Math.exp(getIntercept() + getTimeEffect(i) + residues[i]);
            }
            return rates;
        }

        private double getIntercept() {
            if (intercept) return coefficients.getParameterValue(options.getInterceptIndex());
            else return 0.0;
        }

        private double getTimeEffect(int nodeIndex) {
            if (timeEffect) {

                return coefficients.getParameterValue(options.getTimeEffectIndex()) * nodeTimes[nodeIndex];
            } else {
                return 0.0;
            }
        }

        public double[] getChainSecondDerivative() {
            return orderByTime(getRatesByNodes());
        }

        public double[] getChainGradient() {
            return orderByTime(getRatesByNodes());
        }

        @Override
        public void updateRateGradient(double[] gradient) {
            double[] rates = orderByTime(getRatesByNodes());
            assert(rates.length == gradient.length);
            for (int i = 0; i < gradient.length; i++) {
                gradient[i] *= rates[i];
            }
        }

        @Override
        public void updateRateHessian(double[] gradient, double[] hessian, double[] result) {
            double[] chainGradient = getChainGradient();
            double[] chainSecondDerivative = getChainSecondDerivative();

            for (int i = 0; i < gradient.length; i++) {
                result[i] = gradient[i] * chainSecondDerivative[i] + hessian[i] * chainGradient[i] * chainGradient[i];
            }
        }

        private class Options {

            private final int interceptIndex;
            private final int timeEffectIndex;

            Options (boolean hasIntercept, boolean hasTimeEffect) {
                if (hasTimeEffect) {
                    if (hasIntercept) timeEffectIndex = 1;
                    else timeEffectIndex = 0;
                } else {
                    timeEffectIndex = -1;
                }

                if (hasIntercept) {
                    interceptIndex = 0;
                } else {
                    interceptIndex = -1;
                }
            }

            int getTimeEffectIndex() {
                return timeEffectIndex;
            }

            int getInterceptIndex() {
                return interceptIndex;
            }

        }
    }

}
