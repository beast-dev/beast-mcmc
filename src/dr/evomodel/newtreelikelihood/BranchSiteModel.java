package dr.evomodel.newtreelikelihood;

import dr.evomodel.newsubstmodel.EigenDecomposition;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface BranchSiteModel {

    EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex);

    double[] getStateFrequencies(int categoryIndex);
}
