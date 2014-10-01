package dr.app.beauti.components.dollo;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Marc Suchard
 * @version $Id$
 */

public class DolloComponentFactory implements ComponentFactory {

	public static ComponentFactory INSTANCE = new DolloComponentFactory();

    @Override
    public Class getOptionsClass() {
        return DolloComponentOptions.class;
    }

    public ComponentGenerator createGenerator(BeautiOptions beautiOptions) {
        return new DolloComponentGenerator(beautiOptions);
	}

	public ComponentOptions createOptions(BeautiOptions beautiOptions) {
        return new DolloComponentOptions(beautiOptions);
	}

}