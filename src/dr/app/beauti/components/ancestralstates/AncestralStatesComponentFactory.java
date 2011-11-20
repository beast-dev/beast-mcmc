package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentGenerator;
import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
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
            generator = new SequenceErrorModelComponentGenerator(beautiOptions);
        }
        return generator;
    }

    public ComponentOptions getOptions(final BeautiOptions beautiOptions) {
        if (options == null) {
            options = new SequenceErrorModelComponentOptions();
        }
        return options;
    }

    private SequenceErrorModelComponentGenerator generator = null;
    private SequenceErrorModelComponentOptions options = null;

    public static ComponentFactory INSTANCE = new AncestralStatesComponentFactory();
}
