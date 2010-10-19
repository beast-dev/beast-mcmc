package dr.app.beagle.evomodel.sitemodel;

import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.inference.model.Model;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface BranchSiteModel extends Model {

    EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex);

    double[] getStateFrequencies(int categoryIndex);

    boolean canReturnComplexDiagonalization();
}
