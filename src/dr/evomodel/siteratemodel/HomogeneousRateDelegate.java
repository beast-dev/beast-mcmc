package dr.evomodel.siteratemodel;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

public class HomogeneousRateDelegate extends AbstractModel implements SiteRateDelegate {

    public HomogeneousRateDelegate(String name) {
        super(name);
    }

    @Override
    public int getCategoryCount() {
        return 1;
    }

    @Override
    public void getCategories(double[] categoryRates, double[] categoryProportions) {
        categoryRates[0] = 1.0;
        categoryProportions[0] = 1.0;
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // do nothing
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    }

    protected void acceptState() {
    } // no additional state needs accepting

}
