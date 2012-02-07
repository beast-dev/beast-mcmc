package dr.app.beauti.components.sequenceerror;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.ComponentOptions;
import dr.app.beauti.options.BeautiOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentFactory implements ComponentFactory {

    private SequenceErrorModelComponentFactory() {
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

    public static ComponentFactory INSTANCE = new SequenceErrorModelComponentFactory();
}
