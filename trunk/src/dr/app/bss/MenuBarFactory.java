package dr.app.bss;

import jam.framework.DefaultMenuBarFactory;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class MenuBarFactory extends DefaultMenuBarFactory {

	public MenuBarFactory() {
		
		  registerMenuFactory(new DefaultFileMenuFactory());
		  registerMenuFactory(new DefaultEditMenuFactory());
		  registerMenuFactory(new DefaultHelpMenuFactory());
		  
	}//END: Constructor

}// END: class
