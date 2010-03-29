package dr.app.tracer.application;

import org.virion.jam.framework.*;
import org.virion.jam.mac.*;

import dr.app.util.OSType;

public class TracerMenuBarFactory extends DefaultMenuBarFactory {

    public TracerMenuBarFactory() {
        if (OSType.isMac()) {
            registerMenuFactory(new TracerMacFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
	        registerMenuFactory(new AnalysisMenuFactory());
            registerMenuFactory(new MacWindowMenuFactory());
            registerMenuFactory(new MacHelpMenuFactory());
        } else {
            registerMenuFactory(new TracerDefaultFileMenuFactory());
            registerMenuFactory(new DefaultEditMenuFactory());
	        registerMenuFactory(new AnalysisMenuFactory());
            registerMenuFactory(new DefaultHelpMenuFactory());
        }

    }
}

