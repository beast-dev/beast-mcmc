package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AncestralStatesComponentFactory implements ComponentFactory {

    private AncestralStatesComponentFactory() {
        // singleton pattern - private constructor
    }

    public ComponentGenerator getGenerator(final BeautiOptions beautiOptions) {
        if (generator == null) {
            generator = new AncestralStatesComponentGenerator(beautiOptions);
        }
        return generator;
    }

    public ComponentOptions getOptions(final BeautiOptions beautiOptions) {
        if (options == null) {
            options = new AncestralStatesComponentOptions();
        }
        return options;
    }

    private AncestralStatesComponentGenerator generator = null;
    private AncestralStatesComponentOptions options = null;

    public static ComponentFactory INSTANCE = new AncestralStatesComponentFactory();
}
