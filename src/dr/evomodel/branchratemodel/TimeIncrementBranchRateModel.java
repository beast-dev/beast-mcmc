/*
 * TimeIncrementBranchRateModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 * @author Marc A. Suchard
 * @author Andy Magee
 */

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.List;

public class TimeIncrementBranchRateModel extends AbstractBranchRateModel implements DifferentiableBranchRates, Citable {

    private static final String PARSER_NAME = "timeIncrementBranchRateModel";
    private static final String INCREMENT = "increment";

    private final AbstractBranchRateModel baseModel;
    private final DifferentiableBranchRates baseRate;
    private final int taxonNumber;
    private final double increment;

    public TimeIncrementBranchRateModel(String name,
                                        AbstractBranchRateModel baseModel,
                                        DifferentiableBranchRates baseRate,
                                        Taxon taxon, double increment
    ) {
        super(name);
        this.baseModel = baseModel;
        this.baseRate = baseRate;
        this.taxonNumber = 0; // TODO Fix
        this.increment = increment;

        addModel(baseModel);
    }

    @Override
    public double getUntransformedBranchRate(Tree tree, NodeRef node) {
        return baseRate.getUntransformedBranchRate(tree, node);
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
        return baseRate.getRateParameter();
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return baseRate.getParameterIndexFromNode(node);
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        return baseRate.getTransform();
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
    public double getBranchRate(Tree tree, NodeRef node) {
        return 0; // TODO Figure this out
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARSER_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            RandomLocalClockModel rlcModel = (RandomLocalClockModel) xo.getChild(RandomLocalClockModel.class);
            AbstractBranchRateModel branchRates = (AbstractBranchRateModel) xo.getChild(AbstractBranchRateModel.class);
            if (!(branchRates instanceof DifferentiableBranchRates)) {
                throw new XMLParseException("Need a good comment");
            }
            DifferentiableBranchRates differentiableBranchRates = (DifferentiableBranchRates) branchRates;

            Taxon taxon = (Taxon) xo.getChild(Taxon.class);
            double length = xo.getDoubleAttribute(INCREMENT);

            return new TimeIncrementBranchRateModel("name",
                    branchRates, differentiableBranchRates, taxon, length);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(AbstractBranchRateModel.class),
                new ElementRule(Taxon.class),
                AttributeRule.newDoubleRule(INCREMENT),
        };

        public String getParserDescription() {
            return "TODO"; // TODO
        }

        public Class getReturnType() {
            return TimeIncrementBranchRateModel.class;
        }
    };

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged(object, index);
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new RuntimeException("Should not get here");
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
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        return null;
    }
}
