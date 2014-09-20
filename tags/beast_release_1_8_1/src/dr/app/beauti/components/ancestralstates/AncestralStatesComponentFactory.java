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

    @Override
    public Class getOptionsClass() {
        return AncestralStatesComponentOptions.class;
    }

    @Override
    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new AncestralStatesComponentGenerator(beautiOptions);
    }

    @Override
    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new AncestralStatesComponentOptions();
    }

    public static ComponentFactory INSTANCE = new AncestralStatesComponentFactory();
}
