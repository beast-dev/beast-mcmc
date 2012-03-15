package dr.app.beauti.components.hpm;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public class HierarchicalModelComponentFactory implements ComponentFactory {

    private HierarchicalModelComponentFactory() {
        // singleton pattern - private constructor
    }

    public ComponentGenerator getGenerator(final BeautiOptions beautiOptions) {
        if (generator == null) {
            generator = new HierarchicalModelComponentGenerator(beautiOptions);
        }
        return generator;
    }

    public ComponentOptions getOptions(final BeautiOptions beautiOptions) {
        if (options == null) {
            options = new HierarchicalModelComponentOptions(beautiOptions);
        }
        return options;
    }

    private HierarchicalModelComponentGenerator generator = null;
    private HierarchicalModelComponentOptions options = null;

    public static ComponentFactory INSTANCE = new HierarchicalModelComponentFactory();
}