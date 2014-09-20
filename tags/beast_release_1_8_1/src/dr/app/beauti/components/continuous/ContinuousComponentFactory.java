package dr.app.beauti.components.continuous;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */

public class ContinuousComponentFactory implements ComponentFactory {

	public static ComponentFactory INSTANCE = new ContinuousComponentFactory();

    @Override
    public Class getOptionsClass() {
        return ContinuousComponentOptions.class;
    }

    public ComponentGenerator createGenerator(BeautiOptions beautiOptions) {
		return new ContinuousComponentGenerator(beautiOptions);
	}

	public ComponentOptions createOptions(BeautiOptions beautiOptions) {
        return new ContinuousComponentOptions(beautiOptions);
	}

}