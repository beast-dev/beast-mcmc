package dr.app.beagle.evomodel.sitemodel;

import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Model;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public interface BranchSubstitutionModel extends Model {

    EigenDecomposition getEigenDecomposition(int modelIndex, int categoryIndex);

    double[] getStateFrequencies(int categoryIndex);

    public int getBranchIndex(final Tree tree, final NodeRef node);

    public int getEigenCount();

    boolean canReturnComplexDiagonalization();
}
