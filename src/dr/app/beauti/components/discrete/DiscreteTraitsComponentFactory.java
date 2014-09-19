package dr.app.beauti.components.discrete;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.ComponentOptions;
import dr.app.beauti.options.BeautiOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTraitsComponentFactory implements ComponentFactory {

    private DiscreteTraitsComponentFactory() {
        // singleton pattern - private constructor
    }

    @Override
    public Class getOptionsClass() {
        return DiscreteTraitsComponentOptions.class;
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new DiscreteTraitsComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new DiscreteTraitsComponentOptions(beautiOptions);
    }

    public static ComponentFactory INSTANCE = new DiscreteTraitsComponentFactory();
}