/*
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

package dr.app.beauti.options;

import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.types.StartingTreeType;
import dr.app.beauti.types.TreePriorType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.Tree;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionTreeModel extends PartitionOptions {

    private PartitionTreePrior treePrior;

    private StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    private Tree userStartingTree = null;

    private boolean isNewick = true;
    private boolean fixedTree = false;
//    private double initialRootHeight = 1.0;

    //TODO if use EBSP and *BEAST, validate Ploidy of every PD is same for each tree that the PD(s) belongs to
    // BeastGenerator.checkOptions()
    private PloidyType ploidyType = PloidyType.AUTOSOMAL_NUCLEAR;

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
        fixedTree = source.fixedTree;
//        initialRootHeight = source.initialRootHeight;
        ploidyType = source.ploidyType;
    }

    protected void initModelParametersAndOpererators() {

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameterTree(this, "treeModel.rootHeight", "root height of the tree", true, 1.0);

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

        //=============== microsat ======================
        createParameter("treeModel.microsatellite.internalNodesParameter", "Microsatellite sampler tree internal node parameter");
        createOperator("microsatInternalNodesParameter", "Microsat tree internal node",
                "Random integer walk on microsatellite sampler tree internal node parameter",
                "treeModel.microsatellite.internalNodesParameter", OperatorType.RANDOM_WALK_INT, 1.0, branchWeights);
    }

    /**
     * return a list of parameters that are required
     *
     * @param parameters the parameter list
     */
    public void selectParameters(List<Parameter> parameters) {
        setAvgRootAndRate();

        getParameter("tree");
        getParameter("treeModel.internalNodeHeights");
        getParameter("treeModel.allInternalNodeHeights");

        Parameter rootHeightParameter = getParameter("treeModel.rootHeight");
        if (rootHeightParameter.priorType == PriorType.NONE_TREE_PRIOR || !rootHeightParameter.isPriorEdited()) {
            rootHeightParameter.initial = getInitialRootHeight();
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

        if (getDataType().getType() == DataType.MICRO_SAT) {
            getParameter("treeModel.microsatellite.internalNodesParameter");
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param operators the operator list
     */
    public void selectOperators(List<Operator> operators) {
        setAvgRootAndRate();

        // if not a fixed tree then sample tree space
        if (!fixedTree) {
            Operator subtreeSlideOp = getOperator("subtreeSlide");
            if (!subtreeSlideOp.tuningEdited) {
                subtreeSlideOp.tuning = getInitialRootHeight() / 10.0;
            }

            operators.add(subtreeSlideOp);
            operators.add(getOperator("narrowExchange"));
            operators.add(getOperator("wideExchange"));
            operators.add(getOperator("wilsonBalding"));
        }

        operators.add(getOperator("treeModel.rootHeight"));
        operators.add(getOperator("uniformHeights"));

        if (getDataType().getType() == DataType.MICRO_SAT) {
            operators.add(getOperator("microsatInternalNodesParameter"));
        }
    }

    /////////////////////////////////////////////////////////////

    public PartitionTreePrior getPartitionTreePrior() {
        return treePrior;
    }

    public void setPartitionTreePrior(PartitionTreePrior treePrior) {
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

    public void setNewick(boolean newick) {
        isNewick = newick;
    }

    public void setPloidyType(PloidyType ploidyType) {
        this.ploidyType = ploidyType;
    }

    public PloidyType getPloidyType() {
        return ploidyType;
    }

    public double getInitialRootHeight() {
        return getAvgRootAndRate()[0];
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
}
