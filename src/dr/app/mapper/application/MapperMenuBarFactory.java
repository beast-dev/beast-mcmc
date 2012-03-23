package dr.app.mapper.application;

import jam.framework.*;
import jam.mac.*;

import dr.app.util.OSType;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MapperMenuBarFactory extends DefaultMenuBarFactory {

    public MapperMenuBarFactory() {
        if (OSType.isMac()) {
            registerMenuFactory(new MapperMacFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
            registerMenuFactory(new MacWindowMenuFactory());
            registerMenuFactory(new MacHelpMenuFactory());
        } else {
            registerMenuFactory(new MapperDefaultFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
            registerMenuFactory(new DefaultHelpMenuFactory());
        }

    }
}

