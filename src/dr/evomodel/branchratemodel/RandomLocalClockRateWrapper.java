/*
 * BranchRates.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * @author Alexander Fisher
 */

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.function.DoubleBinaryOperator;

public class RandomLocalClockRateWrapper implements DifferentiableBranchRates {

    public static final String RLC_RATES_WRAPPER = "rlcRatesWrapper";

    private RandomLocalClockModel rlcModel;

    public RandomLocalClockRateWrapper(RandomLocalClockModel rlcModel) {

        this.rlcModel = rlcModel;
    }

    @Override
    public double getUntransformedBranchRate(Tree tree, NodeRef node) {
        return rlcModel.getUnscaledBranchRate(tree, node);
    }


    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getRateParameter() {
        return null;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return 0;
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        return null;
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");

    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");

    }

    @Override
    public double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        return 0;
    }

    @Override
    public void forEachOverRates(NodeRateMap map) {

    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        return 0;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return RLC_RATES_WRAPPER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            RandomLocalClockModel rlcModel = (RandomLocalClockModel) xo.getChild(RandomLocalClockModel.class);

            RandomLocalClockRateWrapper rlcRateWrapper = new RandomLocalClockRateWrapper(rlcModel);
            return rlcRateWrapper;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RandomLocalClockModel.class),
        };

        public String getParserDescription() {
            return "DifferentiableBranchRates object that wraps around randomLocalClock rates (product of rates and indicators)";
        }

        public Class getReturnType() {
            return RandomLocalClockRateWrapper.class;
        }
    };

}
