package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class LinkedParametersComponentFactory implements ComponentFactory {

    private LinkedParametersComponentFactory() {
        // singleton pattern - private constructor
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new LinkedParametersComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new LinkedParametersComponentOptions(beautiOptions);
    }

    public static ComponentFactory INSTANCE = new LinkedParametersComponentFactory();
}