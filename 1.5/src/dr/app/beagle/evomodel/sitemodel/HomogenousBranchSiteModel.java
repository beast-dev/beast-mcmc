package dr.app.beagle.evomodel.sitemodel;

import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class HomogenousBranchSiteModel extends AbstractModel implements BranchSiteModel {
    public HomogenousBranchSiteModel(SubstitutionModel substModel, FrequencyModel frequencyModel) {
        super("HomogenousBranchSiteModel");
        
        this.substModel = substModel;
        addModel(substModel);
        this.frequencyModel = frequencyModel;
        addModel(frequencyModel);
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

    /**
     * Homogenous model - returns if substitution model can return complex diagonalization
     * @return
     */
    public boolean canReturnComplexDiagonalization() {
        return substModel.canReturnComplexDiagonalization();
    }

    private final SubstitutionModel substModel;
    private final FrequencyModel frequencyModel;

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }
}
