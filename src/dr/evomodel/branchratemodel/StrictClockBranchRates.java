/*
 * StrictClockBranchRates.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.function.DoubleBinaryOperator;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: StrictClockBranchRates.java,v 1.3 2006/01/09 17:44:30 rambaut Exp $
 */
public class StrictClockBranchRates extends AbstractBranchRateModel implements DifferentiableBranchRates {

    private final Parameter rateParameter;

    public StrictClockBranchRates(Parameter rateParameter) {

        super(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES);

        this.rateParameter = rateParameter;

        addVariable(rateParameter);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // nothing to do
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        return rateParameter.getParameterValue(0);
    }

    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        return 1.0;
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getRateParameter() { return rateParameter; }

    @Override
    public int getParameterIndexFromNode(NodeRef node) { return node.getNumber(); }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {

        double total = 0.0;
        for (double v : gradient) {
            total += v;
        }

        return new double[] { total };
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient,
                                                    double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }
}
