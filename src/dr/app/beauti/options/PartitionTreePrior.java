package dr.app.beauti.options;

import dr.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class PartitionTreePrior extends ModelOptions {

    // Instance variables

    private final BeautiOptions options;        
    private String name;
    
    private PartitionTreeModel treeModel;
    
	private TreePrior nodeHeightPrior = TreePrior.CONSTANT;
    private int parameterization = GROWTH_RATE;
    private int skylineGroupCount = 10;
    private int skylineModel = CONSTANT_SKYLINE;
    private int skyrideSmoothing = SKYRIDE_TIME_AWARE_SMOOTHING;
    // AR - this seems to be set to taxonCount - 1 so we don't need to
    // have a settable variable...
    // public int skyrideIntervalCount = 1;
//    public String extendedSkylineModel = VariableDemographicModel.LINEAR;
    private boolean multiLoci = false;
    private double birthDeathSamplingProportion = 1.0;
    private boolean fixedTree = false;
    
    
    public PartitionTreePrior(BeautiOptions options, PartitionTreeModel treeModel) {
    	this.options = options;
		this.name = treeModel.getName();
		this.treeModel = treeModel;
    }

    /**
     * A copy constructor
     *
     * @param options the beauti options
     * @param name    the name of the new model
     * @param source  the source model
     */
    public PartitionTreePrior(BeautiOptions options, String name, PartitionTreePrior source) {
    	this.options = options;
		this.name = name;
		this.treeModel = source.treeModel;

		this.nodeHeightPrior = source.nodeHeightPrior;
		this.parameterization = source.parameterization;
		this.skylineGroupCount = source.skylineGroupCount;
		this.skylineModel = source.skylineModel;
		this.skyrideSmoothing = source.skyrideSmoothing;
		this.multiLoci = source.multiLoci;
		this.birthDeathSamplingProportion = source.birthDeathSamplingProportion;
		this.fixedTree = source.fixedTree;
        
    }

//    public PartitionTreePrior(BeautiOptions options, String name) {
//        this.options = options;
//        this.name = name;
//    }    
    

    public Parameter getParameter(String name) {

        if (name.startsWith(getName())) {
            name = name.substring(getName().length() + 1);
        }
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

        operator.setPrefix(getName());

        return operator;
    }

 
    public String getPrefix() {
        String prefix = "";
        if (options.getActivePartitionTreePriors().size() > 1 ) {//|| options.isSpeciesAnalysis()
            // There is more than one active partition model, or doing species analysis
            prefix += getName() + ".";
        }
        return prefix;
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
    
    /////////////////////////////////////////////////////////////////////////

    public PartitionTreeModel getTreeModel() {
		return treeModel;
	}

	public void setTreeModel(PartitionTreeModel treeModel) {
		this.treeModel = treeModel;
	}
   
    public TreePrior getNodeHeightPrior() {
		return nodeHeightPrior;
	}

	public void setNodeHeightPrior(TreePrior nodeHeightPrior) {
		this.nodeHeightPrior = nodeHeightPrior;
	}

	public int getParameterization() {
		return parameterization;
	}

	public void setParameterization(int parameterization) {
		this.parameterization = parameterization;
	}

	public int getSkylineGroupCount() {
		return skylineGroupCount;
	}

	public void setSkylineGroupCount(int skylineGroupCount) {
		this.skylineGroupCount = skylineGroupCount;
	}

	public int getSkylineModel() {
		return skylineModel;
	}

	public void setSkylineModel(int skylineModel) {
		this.skylineModel = skylineModel;
	}

	public int getSkyrideSmoothing() {
		return skyrideSmoothing;
	}

	public void setSkyrideSmoothing(int skyrideSmoothing) {
		this.skyrideSmoothing = skyrideSmoothing;
	}

	public boolean isMultiLoci() {
		return multiLoci;
	}

	public void setMultiLoci(boolean multiLoci) {
		this.multiLoci = multiLoci;
	}

	public double getBirthDeathSamplingProportion() {
		return birthDeathSamplingProportion;
	}

	public void setBirthDeathSamplingProportion(double birthDeathSamplingProportion) {
		this.birthDeathSamplingProportion = birthDeathSamplingProportion;
	}

	public boolean isFixedTree() {
		return fixedTree;
	}

	public void setFixedTree(boolean fixedTree) {
		this.fixedTree = fixedTree;
	}

}