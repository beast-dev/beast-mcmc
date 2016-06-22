/*
 * VariableDemographicModel.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.VariableDemographicModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;

/**
 * @author Joseph Heled
 * @version $Id$
 */
public class VariableDemographicModel extends DemographicModel implements MultiLociTreeSet, Citable {

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
        super(VariableDemographicModelParser.MODEL_NAME);

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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Extended Bayesian Skyline multi-locus coalescent model";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{
                                new Author("J", "Heled"),
                                new Author("AJ", "Drummond"),
                        },
                        "Bayesian inference of population size history from multiple loci",
                        2008,
                        "BMC Evolutionary Biology",
                        8,
                        "289",
                        "10.1186/1471-2148-8-289"
                ));
    }

}
