package dr.app.bss;

import jam.framework.DefaultMenuBarFactory;

public class MenuBarFactory extends DefaultMenuBarFactory {

	public MenuBarFactory() {
		
		  registerMenuFactory(new DefaultFileMenuFactory());
		
	}//END: Constructor
	
//	public void populateMenuBar(JMenuBar menuBar, //
//			DocumentFrame documentFrame, //
//			Application application //
//	) {
//
//		super.populateMenuBar(menuBar, documentFrame);
//
//		documentFrame.getOpenAction().setEnabled(true);
//		
//		documentFrame.getSaveAction().setEnabled(false);
//		documentFrame.getSaveAsAction().setEnabled(false);
//		documentFrame.getCutAction().setEnabled(false);
//		documentFrame.getCopyAction().setEnabled(false);
//		documentFrame.getPasteAction().setEnabled(false);
//		documentFrame.getDeleteAction().setEnabled(false);
//		documentFrame.getSelectAllAction().setEnabled(false);
//		documentFrame.getFindAction().setEnabled(false);
//
//		documentFrame.getZoomWindowAction().setEnabled(false);
//
//	}// END: populateMenuBar

}// END: class
