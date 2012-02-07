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

    public ComponentGenerator getGenerator(final BeautiOptions beautiOptions) {
        if (generator == null) {
            generator = new TipDateSamplingComponentGenerator(beautiOptions);
        }
        return generator;
    }

    public ComponentOptions getOptions(final BeautiOptions beautiOptions) {
        if (options == null) {
            options = new TipDateSamplingComponentOptions(beautiOptions);
        }
        return options;
    }

    private TipDateSamplingComponentGenerator generator = null;
    private TipDateSamplingComponentOptions options = null;

    public static ComponentFactory INSTANCE = new TipDateSamplingComponentFactory();
}