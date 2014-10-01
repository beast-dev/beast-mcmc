package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class LinkedParameterComponentFactory implements ComponentFactory {

    private LinkedParameterComponentFactory() {
        // singleton pattern - private constructor
    }

    @Override
    public Class getOptionsClass() {
        return LinkedParameterComponentOptions.class;
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new LinkedParameterComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new LinkedParameterComponentOptions(beautiOptions);
    }

    public static ComponentFactory INSTANCE = new LinkedParameterComponentFactory();
}