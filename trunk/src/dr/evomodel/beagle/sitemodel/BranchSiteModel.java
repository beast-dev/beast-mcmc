package dr.evomodel.beagle.sitemodel;

import dr.evomodel.beagle.substmodel.EigenDecomposition;
import dr.inference.model.Model;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface BranchSiteModel extends Model {

    EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex);

    double[] getStateFrequencies(int categoryIndex);
}
