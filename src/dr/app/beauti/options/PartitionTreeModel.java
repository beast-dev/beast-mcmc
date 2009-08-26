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

import dr.app.beauti.enumTypes.ClockType;
import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.StartingTreeType;
import dr.app.beauti.enumTypes.PriorType;
import dr.evolution.datatype.PloidyType;
import dr.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionTreeModel extends PartitionOptions {

    // Instance variables

    private final BeautiOptions options;
    private String name;
    private PartitionTreePrior treePrior;
    private List<PartitionData> allPartitionData = new ArrayList<PartitionData>();

    private StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    private Tree userStartingTree = null;
    
    private double initialRootHeight = 1.0;

	private boolean fixedTree = false;

    //TODO if use EBSP and *BEAST, validate Ploidy of every PD is same for each tree that the PD(s) belongs to
    // BeastGenerator.checkOptions()
    private PloidyType ploidyType = PloidyType.AUTOSOMAL_NUCLEAR;


    public PartitionTreeModel(BeautiOptions options, PartitionData partition) {
        this.options = options;
        this.name = partition.getName();
        
        allPartitionData.clear();
        addPartitionData(partition);

        initTreeModelParaAndOpers();
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTreeModel(BeautiOptions options, String name, PartitionTreeModel source) {
        this.options = options;
        this.name = name;

        this.allPartitionData.clear();
        for (PartitionData partition: source.allPartitionData) {
        	this.allPartitionData.add(partition);			
		} 

        this.startingTreeType = source.startingTreeType;
        this.userStartingTree = source.userStartingTree;

        initTreeModelParaAndOpers();
    }

//    public PartitionTreeModel(BeautiOptions options, String name) {
//        this.options = options;
//        this.name = name;
//    }

    private void initTreeModelParaAndOpers() {
        
        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameter("treeModel.rootHeight", "root height of the tree", true, 1.0, 0.0, Double.POSITIVE_INFINITY);

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
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {
    	calculateInitialRootHeightPerTree();
    	
    	getParameter("tree");
    	getParameter("treeModel.internalNodeHeights");
    	getParameter("treeModel.allInternalNodeHeights");    	
    	
    	Parameter rootHeightPara = getParameter("treeModel.rootHeight");
    	rootHeightPara.initial = initialRootHeight; 
    	rootHeightPara.priorEdited = true;
    	if (!options.starBEASTOptions.isSpeciesAnalysis()) {
    		params.add(rootHeightPara);
    	}
    	     
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
    	calculateInitialRootHeightPerTree();

        // if not a fixed tree then sample tree space
        if (!fixedTree) {
        	Operator subtreeSlideOp = getOperator("subtreeSlide");
            if (!subtreeSlideOp.tuningEdited) {
            	subtreeSlideOp.tuning = initialRootHeight / 10.0;
            }
        	
        	ops.add(subtreeSlideOp);
            ops.add(getOperator("narrowExchange"));
            ops.add(getOperator("wideExchange"));
            ops.add(getOperator("wilsonBalding"));
        }
        
        ops.add(getOperator("treeModel.rootHeight"));
        ops.add(getOperator("uniformHeights"));
    }

    /////////////////////////////////////////////////////////////
    
    public boolean containsUncorrelatedRelaxClock() {
        for (PartitionData partition: allPartitionData) {
            PartitionClockModel clockModel = partition.getPartitionClockModel();
            if (clockModel.getClockType() == ClockType.UNCORRELATED_EXPONENTIAL 
                    || clockModel.getClockType() == ClockType.UNCORRELATED_LOGNORMAL) {
                return true;
            }
        }
        return false;
    }

    public List<PartitionData> getAllPartitionData() {
        return allPartitionData;
    }

    public void clearAllPartitionData() {
        this.allPartitionData.clear();
    }

    public void addPartitionData(PartitionData partition) {
        allPartitionData.add(partition);
    }

    public boolean removePartitionData(PartitionData partition) {
        return allPartitionData.remove(partition);
    }

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

    public void setPloidyType(PloidyType ploidyType) {
        this.ploidyType = ploidyType;
    }

    public PloidyType getPloidyType() {
        return ploidyType;
    }
    
    public double getInitialRootHeight() {
		return initialRootHeight;
	}

	public void setInitialRootHeight(double initialRootHeight) {
		this.initialRootHeight = initialRootHeight;
	}
	
	private void calculateInitialRootHeightPerTree() {			
		initialRootHeight = options.clockModelOptions.calculateInitialRootHeightAndRate(allPartitionData) [0];
	}
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
    }

    public Parameter getParameter(String name) {

        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

        return parameter;
    }

    public Operator getOperator(String name) {

        Operator operator = operators.get(name);

        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");

        operator.setPrefix(getPrefix());

        return operator;
    }


    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }


}
