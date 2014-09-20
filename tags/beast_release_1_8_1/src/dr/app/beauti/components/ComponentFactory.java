package dr.app.beauti.components;

import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface ComponentFactory {
    Class getOptionsClass();
    ComponentGenerator createGenerator(BeautiOptions beautiOptions);
    ComponentOptions createOptions(BeautiOptions beautiOptions);
}
