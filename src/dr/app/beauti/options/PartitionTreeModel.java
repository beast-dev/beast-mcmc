package dr.app.beauti.options;

import java.util.List;
import java.util.ArrayList;

import dr.evolution.tree.Tree;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionTreeModel extends ModelOptions {
	
	// Instance variables

    private final BeautiOptions options;
    private String name;
    private PartitionTreePrior treePrior;
    
    private StartingTreeType startingTreeType = StartingTreeType.RANDOM;
    private Tree userStartingTree = null;    
    
	public PartitionTreeModel(BeautiOptions options, PartitionData partition) {
		this.options = options;
		this.name = partition.getName();
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
		
		this.startingTreeType = source.startingTreeType;
		this.userStartingTree = source.userStartingTree;         
    }

//    public PartitionTreeModel(BeautiOptions options, String name) {
//        this.options = options;
//        this.name = name;
//    }


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

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public List<Operator> getOperators() {
        List<Operator> operators = new ArrayList<Operator>();

        return operators;
    }

    public String toString() {
        return getName();
    }

    /**
     * @param includeRelativeRates true if relative rate parameters should be added
     * @return a list of parameters that are required
     */
    List<Parameter> getParameters(boolean includeRelativeRates) {

        List<Parameter> params = new ArrayList<Parameter>();


        return params;
    }

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
        if (options.getActivePartitionTreeModels().size() > 1) { //|| options.isSpeciesAnalysis()
            // There is more than one active partition model
            prefix += getName() + ".";
        }
        return prefix;
    }


}
