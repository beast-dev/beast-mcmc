package dr.evomodel.newtreelikelihood;

import dr.evomodel.newsubstmodel.EigenDecomposition;
import dr.evomodel.newsubstmodel.SubstitutionModel;
import dr.evomodel.newsubstmodel.FrequencyModel;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class HomogenousBranchSiteModel implements BranchSiteModel {
    public HomogenousBranchSiteModel(SubstitutionModel substModel, FrequencyModel frequencyModel) {
        this.substModel = substModel;
        this.frequencyModel = frequencyModel;
    }

    /**
     * Homogenous model - returns the same substitution model for all branches/categories
     * @param branchIndex
     * @param categoryIndex
     * @return
     */
    public EigenDecomposition getEigenDecomposition(int branchIndex, int categoryIndex) {
        return substModel.getEigenDecomposition();
    }

    /**
     * Homogenous model - returns the same frequency model for all categories
     * @param categoryIndex
     * @return
     */
    public double[] getStateFrequencies(int categoryIndex) {
        return frequencyModel.getFrequencies();
    }

    private final SubstitutionModel substModel;
    private final FrequencyModel frequencyModel;

}
