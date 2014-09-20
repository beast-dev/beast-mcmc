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

    @Override
    public Class getOptionsClass() {
        return HierarchicalModelComponentOptions.class;
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new HierarchicalModelComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new HierarchicalModelComponentOptions(beautiOptions);
    }

    public static ComponentFactory INSTANCE = new HierarchicalModelComponentFactory();
}