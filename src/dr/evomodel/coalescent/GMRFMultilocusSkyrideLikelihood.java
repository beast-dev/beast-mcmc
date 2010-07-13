package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Erik W. Bloomquist
 */

public class GMRFMultilocusSkyrideLikelihood extends GMRFSkyrideLikelihood implements MultiLociTreeSet {

    	public GMRFMultilocusSkyrideLikelihood(List<Tree> trees,
                                               Parameter popParameter,
                                               Parameter groupParameter,
                                               Parameter precParameter,
	                                           Parameter lambda,
                                               Parameter beta,
                                               MatrixParameter dMatrix,
	                                           boolean timeAwareSmoothing) {
            super(trees, popParameter, groupParameter, precParameter, lambda, beta, dMatrix, timeAwareSmoothing);
        }

    protected void setTree(List<Tree> treeList) {
        treesSet = this;
        this.treeList = treeList;
        intervalsList = new ArrayList<TreeIntervals>();
        for (Tree tree : treeList) {
            intervalsList.add(new TreeIntervals(tree));
            if (tree instanceof TreeModel) {
                addModel((TreeModel) tree);
            }
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model instanceof TreeModel) {
            TreeModel treeModel = (TreeModel) model;
            int tn = treeList.indexOf(treeModel);
            if (tn >= 0) {
                intervalsList.get(tn).setIntervalsUnknown();
            }
        }
        
        super.handleModelChangedEvent(model, object, index);
    }

    public void initializationReport() {
		System.out.println("Creating a GMRF smoothed skyride model for multiple loci:");
		System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
		System.out.println("\tIf you publish results using this model, please reference: ");
        System.out.println("\t\tMinin, Bloomquist and Suchard (2008) Molecular Biology and Evolution, 25, 1459-1471, and");
        System.out.println("\t\tSuchard, Drummond and Lemey (in preparation).");
	}

    public void wrapSetupIntervals() {
        // Do nothing
    }

    protected void setupSufficientStatistics() {
        // TODO Much to do here
    }

    private List<Tree> treeList;
    private List<TreeIntervals> intervalsList;

    public int nLoci() {
        return treeList.size();
    }

    public Tree getTree(int nt) {
        return treeList.get(nt);
    }

    public TreeIntervals getTreeIntervals(int nt) {
        return intervalsList.get(nt);
    }

    public double getPopulationFactor(int nt) {
        return 1.0;
    }

    public void storeTheState() {
        for (TreeIntervals intervals : intervalsList) {
            intervals.storeState();
        }
    }

    public void restoreTheState() {
        for (TreeIntervals intervals : intervalsList) {
            intervals.restoreState();
        }
    }
}

