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

    @Override
    public Class getOptionsClass() {
        return SequenceErrorModelComponentOptions.class;
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new SequenceErrorModelComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new SequenceErrorModelComponentOptions();
    }

    public static ComponentFactory INSTANCE = new SequenceErrorModelComponentFactory();
}
