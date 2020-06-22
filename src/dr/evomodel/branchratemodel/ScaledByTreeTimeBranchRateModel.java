/*
 * RandomLocalClockModel.java
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
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.ScaledByTreeTimeBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Alexei J. Drummond
 * @author Alexander Fisher
 */

public class ScaledByTreeTimeBranchRateModel extends AbstractBranchRateModel implements Citable {

    private final TreeModel treeModel;
    private final BranchRateModel branchRateModel;
    private final Parameter meanRateParameter;

    private boolean scaleFactorKnown;
    private boolean storedScaleFactorKnown;

    private double scaleFactor;
    private double storedScaleFactor;

    public ScaledByTreeTimeBranchRateModel(TreeModel treeModel,
                                           BranchRateModel branchRateModel,
                                           Parameter meanRateParameter) {

        super(ScaledByTreeTimeBranchRateModelParser.TREE_TIME_BRANCH_RATES);

        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.meanRateParameter = meanRateParameter;

        addModel(treeModel);
        addModel(branchRateModel);

        if (meanRateParameter != null) {
            addVariable(meanRateParameter);
        }

        scaleFactorKnown = false;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        scaleFactorKnown = false;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        scaleFactorKnown = false;
        fireModelChanged();
    }

    protected void storeState() {
        storedScaleFactor = scaleFactor;
        storedScaleFactorKnown = scaleFactorKnown;
    }

    protected void restoreState() {
        scaleFactor = storedScaleFactor;
        scaleFactorKnown = storedScaleFactorKnown;
    }

    protected void acceptState() { }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == treeModel;

        if (!scaleFactorKnown) {
            scaleFactor = calculateScaleFactor();
        }

        return scaleFactor * branchRateModel.getBranchRate(tree, node);
    }

    private double calculateScaleFactor() {

        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {

                double branchTime = treeModel.getBranchLength(node);

                double branchLength = branchTime * branchRateModel.getBranchRate(treeModel, node);

                timeTotal += branchTime;
                branchTotal += branchLength;
            }
        }

        double scaleFactor = timeTotal / branchTotal;

        if (meanRateParameter != null) {
            scaleFactor *= meanRateParameter.getParameterValue(0);
        }

        return scaleFactor;
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

        description += " with scaling-by-tree-time";
        return description;
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> list = 
                (branchRateModel instanceof Citable) ?
                        new ArrayList<Citation>(((Citable) branchRateModel).getCitations()) :
                        new ArrayList<Citation>();
        list.add(RandomLocalClockModel.CITATION);
        return list;
    }
}