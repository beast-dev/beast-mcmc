/*
 * PartitionTreeModel.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.*;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.Tree;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionTreeModel extends PartitionOptions {

    private static final long serialVersionUID = 4829401415152235625L;

    private PartitionTreePrior treePrior;

    private StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    private Tree userStartingTree = null;

    private boolean isNewick = true;

    //TODO if use EBSP and *BEAST, validate Ploidy of every PD is same for each tree that the PD(s) belongs to
    // BeastGenerator.checkOptions()
    private PloidyType ploidyType = PloidyType.AUTOSOMAL_NUCLEAR;

    private boolean hasTipCalibrations = false;
    private boolean hasNodeCalibrations = false;


    public PartitionTreeModel(BeautiOptions options, AbstractPartitionData partition) {
        super(options, partition.getName());
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTreeModel(BeautiOptions options, String name, PartitionTreeModel source) {
        super(options, name);

        treePrior = source.treePrior;
        startingTreeType = source.startingTreeType;
        userStartingTree = source.userStartingTree;

        isNewick = source.isNewick;
//        initialRootHeight = source.initialRootHeight;
        ploidyType = source.ploidyType;
    }

    protected void initModelParametersAndOpererators() {

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameterTree(this, "treeModel.rootHeight", "root height of the tree", true);

        //TODO treeBitMove should move to PartitionClockModelTreeModelLink, after Alexei finish
        createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree",
                OperatorType.TREE_BIT_MOVE, -1.0, treeWeights);

        createScaleOperator("treeModel.rootHeight", demoTuning, demoWeights);
        createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally",
                "treeModel.internalNodeHeights", OperatorType.UNIFORM, -1, branchWeights);

        createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree",
                OperatorType.SUBTREE_SLIDE, 1.0, treeWeights);
        createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree",
                OperatorType.NARROW_EXCHANGE, -1, treeWeights);
        createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree",
                OperatorType.WIDE_EXCHANGE, -1, demoWeights);
        createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree",
                OperatorType.WILSON_BALDING, -1, demoWeights);

        createOperator("subtreeLeap", "Tree", "Performs the subtree-leap rearrangement of the tree", "tree",
                OperatorType.SUBTREE_LEAP, 1.0, options.taxonList.getTaxonCount() < treeWeights ? treeWeights : options.taxonList.getTaxonCount());

    }

    /**
     * return a list of parameters that are required
     *
     * @param parameters the parameter list
     */
    public void selectParameters(List<Parameter> parameters) {
//        setAvgRootAndRate();

        getParameter("tree");
        getParameter("treeModel.internalNodeHeights");
        getParameter("treeModel.allInternalNodeHeights");

        Parameter rootHeightParameter = getParameter("treeModel.rootHeight");
        if (rootHeightParameter.priorType == PriorType.NONE_TREE_PRIOR || !rootHeightParameter.isPriorEdited()) {
            rootHeightParameter.setInitial(getInitialRootHeight());
            rootHeightParameter.truncationLower = options.maximumTipHeight;
            rootHeightParameter.uniformLower = options.maximumTipHeight;
            rootHeightParameter.isTruncated = true;
        }

        if (options.useStarBEAST) {
            rootHeightParameter.isCalibratedYule = treePrior.getNodeHeightPrior() == TreePriorType.SPECIES_YULE_CALIBRATION;
        } else {
            rootHeightParameter.isCalibratedYule = treePrior.getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION;
            parameters.add(rootHeightParameter);
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param operators the operator list
     */
    public void selectOperators(List<Operator> operators) {
//        setAvgRootAndRate();

        Operator subtreeSlideOp = getOperator("subtreeSlide");
        if (!subtreeSlideOp.isTuningEdited()) {
            double tuning = 1.0;
            if (!Double.isNaN(getInitialRootHeight()) && !Double.isInfinite(getInitialRootHeight())) {
                tuning = getInitialRootHeight() / 10.0;
            }
            subtreeSlideOp.setTuning(tuning);
        }

        operators.add(subtreeSlideOp);
        operators.add(getOperator("narrowExchange"));
        operators.add(getOperator("wideExchange"));
        operators.add(getOperator("wilsonBalding"));

        operators.add(getOperator("treeModel.rootHeight"));
        operators.add(getOperator("uniformHeights"));

        operators.add(getOperator("subtreeLeap"));

        boolean defaultInUse;
        boolean branchesInUse;
        boolean newMixInUse;

        // if not a fixed tree then sample tree space
        if (options.operatorSetType == OperatorSetType.DEFAULT) {
            defaultInUse = true;
            branchesInUse = true;
            newMixInUse = false;
        } else if (options.operatorSetType == OperatorSetType.NEW_TREE_MIX) {
            defaultInUse = false;
            branchesInUse = false;
            newMixInUse = true;
        } else if (options.operatorSetType == OperatorSetType.FIXED_TREE_TOPOLOGY) {
            defaultInUse = false;
            branchesInUse = true;
            newMixInUse = false;
        } else {
            throw new IllegalArgumentException("Unknown operator set type");
        }

        getOperator("subtreeSlide").setUsed(defaultInUse);
        getOperator("narrowExchange").setUsed(defaultInUse);
        getOperator("wideExchange").setUsed(defaultInUse);
        getOperator("wilsonBalding").setUsed(defaultInUse);

        getOperator("treeModel.rootHeight").setUsed(branchesInUse);
        getOperator("uniformHeights").setUsed(branchesInUse);

        getOperator("subtreeLeap").setUsed(newMixInUse);
    }

    /////////////////////////////////////////////////////////////

    public PartitionTreePrior getPartitionTreePrior() {
        return treePrior;
    }

    public void setPartitionTreePrior(PartitionTreePrior treePrior) {
        options.clearDataPartitionCaches();
        this.treePrior = treePrior;
    }

    public StartingTreeType getStartingTreeType() {
        return startingTreeType;
    }

    public void setStartingTreeType(StartingTreeType startingTreeType) {
        this.startingTreeType = startingTreeType;
    }

    public Tree getUserStartingTree() {
        return userStartingTree;
    }

    public void setUserStartingTree(Tree userStartingTree) {
        this.userStartingTree = userStartingTree;
    }

    public boolean isNewick() {
        return isNewick;
    }

    public void setNewick(boolean isNewick) {
        this.isNewick = isNewick;
    }

    public void setTipCalibrations(boolean hasTipCalibrations) {
        this.hasTipCalibrations = hasTipCalibrations;
    }

    public void setNodeCalibrations(boolean hasNodeCalibrations) {
        this.hasNodeCalibrations = hasNodeCalibrations;
    }


    public void setPloidyType(PloidyType ploidyType) {
        this.ploidyType = ploidyType;
    }

    public PloidyType getPloidyType() {
        return ploidyType;
    }

    public double getInitialRootHeight() {
        return Double.NaN;
//        return getAvgRootAndRate()[0];
    }

//    public void setInitialRootHeight(double initialRootHeight) {
//        this.initialRootHeight = initialRootHeight;
//    }

//    private void calculateInitialRootHeightPerTree() {
//		initialRootHeight = options.clockModelOptions
//                .calculateInitialRootHeightAndRate(options.getDataPartitions(this)) [0];
//    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }

    public int getDimension() { // n-1
        return options.getTaxonCount(options.getDataPartitions(this)) - 1;
    }

    public int getTaxonCount() {
        return options.getTaxonCount(options.getDataPartitions(this));
    }
}
