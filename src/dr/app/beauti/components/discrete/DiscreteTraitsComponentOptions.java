package dr.app.beauti.components.discrete;

import dr.app.beauti.options.*;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTraitsComponentOptions implements ComponentOptions {

    private final BeautiOptions options;

    public DiscreteTraitsComponentOptions(final BeautiOptions options) {
        this.options = options;
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createParameter("trait.frequencies", "the frequencies of each state");
        modelOptions.createScaleOperator("trait.frequencies", 0.75, 1.0);
        modelOptions.createParameter("trait.rates", "the instantaneous transition rates between states");
        modelOptions.createScaleOperator("trait.rates", 0.75, 1.0);
        modelOptions.createParameter("trait.indicators", "a vector of bits indicating non-zero rates");
        modelOptions.createScaleOperator("trait.indicators", 0.75, 1.0);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        params.add(modelOptions.getParameter("trait.frequencies"));
        params.add(modelOptions.getParameter("trait.rates"));
        if (usingBSSVS) {
            params.add(modelOptions.getParameter("trait.indicators"));
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        ops.add(modelOptions.getOperator("trait.frequencies"));
        ops.add(modelOptions.getOperator("trait.rates"));
        if (usingBSSVS) {
            ops.add(modelOptions.getOperator("trait.indicators"));
        }
    }


    public String discreteTraitName = "location";
    public boolean usingBSSVS = true;
}