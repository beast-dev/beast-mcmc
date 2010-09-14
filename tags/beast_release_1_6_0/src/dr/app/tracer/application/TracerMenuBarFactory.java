package dr.app.tracer.application;

import jam.framework.*;
import jam.mac.*;

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

