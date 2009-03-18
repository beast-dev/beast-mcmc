package dr.evomodel.beagle.sitemodel;

import dr.evomodel.beagle.substmodel.EigenDecomposition;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface BranchSiteModel {

    EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex);

    double[] getStateFrequencies(int categoryIndex);
}
