package dr.app.bss;

//import jam.framework.DefaultEditMenuFactory;
//import jam.framework.DefaultHelpMenuFactory;
import jam.framework.DefaultMenuBarFactory;

public class MenuBarFactory extends DefaultMenuBarFactory {

	public MenuBarFactory() {
		
		  registerMenuFactory(new DefaultFileMenuFactory());
		  //TODO: override Edit & Help menus
		  registerMenuFactory(new DefaultEditMenuFactory());
		  registerMenuFactory(new DefaultHelpMenuFactory());
		  
	}//END: Constructor

}// END: class
