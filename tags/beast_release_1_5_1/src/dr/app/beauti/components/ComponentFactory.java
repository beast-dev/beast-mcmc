package dr.app.beauti.components;

import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface ComponentFactory {
    ComponentGenerator getGenerator(BeautiOptions beautiOptions);
    ComponentOptions getOptions(BeautiOptions beautiOptions);
}
