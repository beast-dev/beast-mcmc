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

    public ComponentGenerator getGenerator(final BeautiOptions beautiOptions) {
        if (generator == null) {
            generator = new DiscreteTraitsComponentGenerator(beautiOptions);
        }
        return generator;
    }

    public ComponentOptions getOptions(final BeautiOptions beautiOptions) {
        if (options == null) {
            options = new DiscreteTraitsComponentOptions(beautiOptions);
        }
        return options;
    }

    private DiscreteTraitsComponentGenerator generator = null;
    private DiscreteTraitsComponentOptions options = null;

    public static ComponentFactory INSTANCE = new DiscreteTraitsComponentFactory();
}