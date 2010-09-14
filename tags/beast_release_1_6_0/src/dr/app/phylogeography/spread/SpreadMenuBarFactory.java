package dr.app.phylogeography.spread;

import jam.framework.DefaultEditMenuFactory;
import jam.framework.DefaultHelpMenuFactory;
import jam.framework.DefaultMenuBarFactory;
import jam.mac.MacHelpMenuFactory;
import jam.mac.MacWindowMenuFactory;
import jam.mac.Utils;
import dr.app.phylogeography.spread.layerspanel.DefaultFileMenuFactory;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SpreadMenuBarFactory extends DefaultMenuBarFactory {

    public SpreadMenuBarFactory() {
        if (Utils.isMacOSX()) {
            registerMenuFactory(new MacFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
            registerMenuFactory(new MacWindowMenuFactory());
            registerMenuFactory(new MacHelpMenuFactory());
        } else {
            registerMenuFactory(new DefaultFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
            registerMenuFactory(new DefaultHelpMenuFactory());
        }
    }

}