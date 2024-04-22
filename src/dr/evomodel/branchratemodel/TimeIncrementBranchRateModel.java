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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleBinaryOperator;

public class TimeIncrementBranchRateModel extends AbstractBranchRateModel implements DifferentiableBranchRates, Citable {

    private static final String PARSER_NAME = "timeIncrementBranchRateModel";
    private static final String INCREMENT = "increment";

    private final TreeModel treeModel;
    private final BranchRateModel branchRateModel;
    private final DifferentiableBranchRates differentiableBranchRateModel;
    private final NodeRef tip;
    private final Parameter offset;

    public TimeIncrementBranchRateModel(String name,
                                        TreeModel treeModel,
                                        BranchRateModel branchRateModel,
                                        Taxon taxon, Parameter offset) {
        super(name);

        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.differentiableBranchRateModel = (branchRateModel instanceof DifferentiableBranchRates) ?
                (DifferentiableBranchRates) branchRateModel : null;

        this.tip = treeModel.getNode(treeModel.getTaxonIndex(taxon));

        if (tip == null) {
            throw new IllegalArgumentException("Unable to find tip node for taxon");
        }

        this.offset = offset;

        addModel(treeModel); // TODO Can probably remove?
        addModel(branchRateModel);
        addVariable(offset);
    }

    @Override
    public void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
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
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        checkDifferentiability();
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        checkDifferentiability();
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getRateParameter() {
        checkDifferentiability();
        return differentiableBranchRateModel.getRateParameter();
    }

    public Tree getTree() {
        return treeModel;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        checkDifferentiability();
        return differentiableBranchRateModel.getParameterIndexFromNode(node);
    }

    private void checkDifferentiability() {
        if (differentiableBranchRateModel == null) {
            throw new RuntimeException("Non-differentiable base BranchRateModel");
        }
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value,
                                                    int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        checkDifferentiability();
        return differentiableBranchRateModel.mapReduceOverRates(map, reduce, initial);
    }

    @Override
    public void forEachOverRates(NodeRateMap map) {
        checkDifferentiability();
        differentiableBranchRateModel.forEachOverRates(map);
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {
        double rate = branchRateModel.getBranchRate(tree, node);

        if (node == tip) {
            double oldLength = tree.getBranchLength(node);
            double newLength = oldLength + offset.getParameterValue(0);

            rate *= (newLength / oldLength);
        }

        return rate;
    }

    public double getUntransformedBranchRate(Tree tree, NodeRef node) {
        checkDifferentiability();
        throw new RuntimeException("Not yet implemented");
//        return differentiableBranchRateModel.getUntransformedBranchRate(tree, node);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        String description =
                (branchRateModel instanceof Citable) ?
                        ((Citable) branchRateModel).getDescription() :
                        "Unknown clock model";

        description += " with lost time for taxon " + tip;
        return description;
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> list =
                (branchRateModel instanceof Citable) ?
                        new ArrayList<>(((Citable) branchRateModel).getCitations()) :
                        new ArrayList<>();
        // TODO
        return list;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARSER_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);;
            BranchRateModel branchRates = (BranchRateModel) xo.getChild(AbstractBranchRateModel.class);
            Taxon taxon = (Taxon) xo.getChild(Taxon.class);
            Parameter offset = (Parameter) xo.getChild(Parameter.class);

            return new TimeIncrementBranchRateModel("name", tree, branchRates, taxon, offset);
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
                new ElementRule(Parameter.class),
        };

        public String getParserDescription() {
            return "TODO"; // TODO
        }

        public Class getReturnType() {
            return TimeIncrementBranchRateModel.class;
        }
    };
}
