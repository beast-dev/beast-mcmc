/*
 * VariableDemographicModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Joseph Heled
 * @version $Id$
 */
public class VariableDemographicModel extends DemographicModel implements MultiLociTreeSet {
    static final String MODEL_NAME = "variableDemographic";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String INDICATOR_PARAMETER = "indicators";
    public static final String POPULATION_TREES = "trees";
    private static final String PLOIDY = SpeciesBindings.PLOIDY;
    public static String POP_TREE = "ptree";

    public static final String LOG_SPACE = "logUnits";
    public static final String USE_MIDPOINTS = "useMidpoints";

    public static final String TYPE = "type";
    //public static final String STEPWISE = "stepwise";
    //public static final String LINEAR = "linear";
    //public static final String EXPONENTIAL = "exponential";

    public static final String demoElementName = "demographic";

    private final Parameter popSizeParameter;
    private final Parameter indicatorParameter;

    public Type getType() {
        return type;
    }

    private final Type type;
    private final boolean logSpace;
    private final boolean mid;
    private final TreeModel[] trees;
    private VDdemographicFunction demoFunction = null;
    private VDdemographicFunction savedDemoFunction = null;
    private final double[] populationFactors;

    public Parameter getPopulationValues() {
        return popSizeParameter;
    }

    public enum Type {
        LINEAR("linear"),
        EXPONENTIAL("exponential"),
        STEPWISE("stepwise");

        Type(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        String name;
    }

    public VariableDemographicModel(TreeModel[] trees, double[] popFactors,
                                    Parameter popSizeParameter, Parameter indicatorParameter,
                                    Type type, boolean logSpace, boolean mid) {
        super(MODEL_NAME);

        this.popSizeParameter = popSizeParameter;
        this.indicatorParameter = indicatorParameter;

        this.populationFactors = popFactors;

        int events = 0;
        for (Tree t : trees) {
            // number of coalescent events
            events += t.getExternalNodeCount() - 1;
            // we will have to handle this I guess
            assert t.getUnits() == trees[0].getUnits();
        }
        // all trees share time 0, need fixing for serial data

        events += type == Type.STEPWISE ? 0 : 1;

        final int popSizes = popSizeParameter.getDimension();
        final int nIndicators = indicatorParameter.getDimension();
        this.type = type;
        this.logSpace = logSpace;
        this.mid = mid;

        if (popSizes != events) {

            System.err.println("INFO: resetting length of parameter " + popSizeParameter.getParameterName() +
                    "(size " + popSizeParameter.getSize() + ") in variable demographic model to " + events);
            popSizeParameter.setDimension(events);
            popSizeParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, -Double.MAX_VALUE, popSizeParameter.getDimension()));
        }

        if (nIndicators != events - 1) {
            System.err.println("INFO: resetting length of parameter " + indicatorParameter.getParameterName() +
                    " in variable demographic model to " + (events - 1));

            indicatorParameter.setDimension(events - 1);
            indicatorParameter.addBounds(new Parameter.DefaultBounds(1, 0, indicatorParameter.getDimension()));
        }

        this.trees = trees;

        for (TreeModel t : trees) {
            addModel(t);
        }

        addVariable(indicatorParameter);
        addVariable(popSizeParameter);
    }

    public int nLoci() {
        return trees.length;
    }

    public Tree getTree(int k) {
        return trees[k];
    }

    public TreeIntervals getTreeIntervals(int nt) {
        return getDemographicFunction().getTreeIntervals(nt);
    }

    public double getPopulationFactor(int nt) {
        return populationFactors[nt];
    }

    public void storeTheState() {
        // as a demographic model store/restore is already taken care of 
    }

    public void restoreTheState() {
        // as a demographic model store/restore is already taken care of
    }

    public VDdemographicFunction getDemographicFunction() {
        if (demoFunction == null) {
            demoFunction = new VDdemographicFunction(trees, type,
                    indicatorParameter.getParameterValues(), popSizeParameter.getParameterValues(), logSpace, mid);
        } else {
            demoFunction.setup(trees, indicatorParameter.getParameterValues(), popSizeParameter.getParameterValues(),
                    logSpace, mid);
        }
        return demoFunction;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // tree has changed
        //System.out.println("model changed: " + model);
        if (demoFunction != null) {
            if (demoFunction == savedDemoFunction) {
                demoFunction = new VDdemographicFunction(demoFunction);
            }
            for (int k = 0; k < trees.length; ++k) {
                if (model == trees[k]) {
                    demoFunction.treeChanged(k);
                    //System.out.println("tree changed: " + k + " " + Arrays.toString(demoFunction.dirtyTrees)
                    //       + " " + demoFunction.dirtyTrees);
                    break;
                }
                assert k + 1 < trees.length;
            }
        }
        super.handleModelChangedEvent(model, object, index);
        fireModelChanged(this);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        //System.out.println("parm changed: " + parameter);
        super.handleVariableChangedEvent(variable, index, type);
        if (demoFunction != null) {
            if (demoFunction == savedDemoFunction) {
                demoFunction = new VDdemographicFunction(demoFunction);
            }
            demoFunction.setDirty();
        }
        fireModelChanged(this);
    }

    protected void storeState() {
        savedDemoFunction = demoFunction;
    }

    protected void restoreState() {
        //System.out.println("restore");
        demoFunction = savedDemoFunction;
        savedDemoFunction = null;
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return VariableDemographicModel.MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(VariableSkylineLikelihood.POPULATION_SIZES);
            Parameter popParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(VariableSkylineLikelihood.INDICATOR_PARAMETER);
            Parameter indicatorParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(POPULATION_TREES);

            final int nc = cxo.getChildCount();
            TreeModel[] treeModels = new TreeModel[nc];
            double[] populationFactor = new double[nc];

            for (int k = 0; k < treeModels.length; ++k) {
                final XMLObject child = (XMLObject) cxo.getChild(k);
                populationFactor[k] = child.hasAttribute(PLOIDY) ? child.getDoubleAttribute(PLOIDY) : 1.0;

                treeModels[k] = (TreeModel) child.getChild(TreeModel.class);
            }

            Type type = Type.STEPWISE;

            if (xo.hasAttribute(TYPE)) {
                final String s = xo.getStringAttribute(TYPE);
                if (s.equalsIgnoreCase(Type.STEPWISE.toString())) {
                    type = Type.STEPWISE;
                } else if (s.equalsIgnoreCase(Type.LINEAR.toString())) {
                    type = Type.LINEAR;
                } else if (s.equalsIgnoreCase(Type.EXPONENTIAL.toString())) {
                    type = Type.EXPONENTIAL;
                } else {
                    throw new XMLParseException("Unknown Bayesian Skyline type: " + s);
                }
            }

            final boolean logSpace = xo.getAttribute(LOG_SPACE, false) || type == Type.EXPONENTIAL;
            final boolean useMid = xo.getAttribute(USE_MIDPOINTS, false);

            Logger.getLogger("dr.evomodel").info("Variable demographic: " + type.toString() + " control points");

            return new VariableDemographicModel(treeModels, populationFactor, popParam, indicatorParam, type,
                    logSpace, useMid);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the population size vector.";
        }

        public Class getReturnType() {
            return DemographicModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(VariableSkylineLikelihood.TYPE, true),
                AttributeRule.newBooleanRule(LOG_SPACE, true),
                AttributeRule.newBooleanRule(USE_MIDPOINTS, true),

                new ElementRule(VariableSkylineLikelihood.POPULATION_SIZES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(VariableSkylineLikelihood.INDICATOR_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(POPULATION_TREES, new XMLSyntaxRule[]{
                        new ElementRule(POP_TREE, new XMLSyntaxRule[]{
                                AttributeRule.newDoubleRule(PLOIDY, true),
                                new ElementRule(TreeModel.class),
                        }, 1, Integer.MAX_VALUE)
                })
        };
    };
}
