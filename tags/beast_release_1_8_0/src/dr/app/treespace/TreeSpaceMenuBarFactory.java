package dr.app.treespace;

import jam.framework.DefaultEditMenuFactory;
import jam.framework.DefaultFileMenuFactory;
import jam.framework.DefaultHelpMenuFactory;
import jam.framework.DefaultMenuBarFactory;
import jam.mac.MacHelpMenuFactory;
import jam.mac.MacWindowMenuFactory;
import jam.mac.Utils;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TreeSpaceMenuBarFactory extends DefaultMenuBarFactory {

    public TreeSpaceMenuBarFactory() {
        if (Utils.isMacOSX()) {
            registerMenuFactory(new MacFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
            registerMenuFactory(new MacWindowMenuFactory());
            registerMenuFactory(new MacHelpMenuFactory());
        } else {
            registerMenuFactory(new DefaultFileMenuFactory(true));
            registerMenuFactory(new DefaultEditMenuFactory());
            registerMenuFactory(new DefaultHelpMenuFactory());
        }
    }

}