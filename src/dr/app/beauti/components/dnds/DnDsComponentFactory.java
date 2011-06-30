package dr.app.beauti.components.dnds;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Filip Bielejec
 * @version $Id$
 */

public class DnDsComponentFactory implements ComponentFactory {

	public static ComponentFactory INSTANCE = new DnDsComponentFactory();

	private DnDsComponentGenerator generator = null;
	private DnDsComponentOptions options = null;

	public ComponentGenerator getGenerator(BeautiOptions beautiOptions) {
		if (generator == null) {
			generator = new DnDsComponentGenerator(beautiOptions);
		}
		return generator;
	}

	public ComponentOptions getOptions(BeautiOptions beautiOptions) {
		if (options == null) {
			options = new DnDsComponentOptions(beautiOptions);
		}
		return options;
	}

}
