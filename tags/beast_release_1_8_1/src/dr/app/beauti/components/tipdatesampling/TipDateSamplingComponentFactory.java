package dr.app.beauti.components.tipdatesampling;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.ComponentOptions;
import dr.app.beauti.options.BeautiOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TipDateSamplingComponentFactory implements ComponentFactory {

    private TipDateSamplingComponentFactory() {
        // singleton pattern - private constructor
    }

    @Override
    public Class getOptionsClass() {
        return TipDateSamplingComponentOptions.class;
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new TipDateSamplingComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new TipDateSamplingComponentOptions(beautiOptions);
    }

    public static ComponentFactory INSTANCE = new TipDateSamplingComponentFactory();
}