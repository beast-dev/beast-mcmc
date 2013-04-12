package dr.app.bss;

import java.awt.event.KeyEvent;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.text.DefaultEditorKit;

import jam.framework.AbstractFrame;
import jam.framework.MenuBarFactory;
import jam.framework.MenuFactory;

public class DefaultEditMenuFactory implements MenuFactory  {

	public DefaultEditMenuFactory() {
	}// END: Constructor
	
	@Override
	public void populateMenu(JMenu menu, AbstractFrame frame) {

		menu.setMnemonic('E');
		
		JMenuItem item;
		
		// Setup Cut
		item = new JMenuItem(new DefaultEditorKit.CutAction());
		item.setText("Cut");
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_X,
				MenuBarFactory.MENU_MASK));
		menu.add(item);
		
		// Setup Copy
		item = new JMenuItem(new DefaultEditorKit.CopyAction());
		item.setText("Copy");
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C,
				MenuBarFactory.MENU_MASK));
		menu.add(item);
		
		// Setup Paste
		item = new JMenuItem(new DefaultEditorKit.PasteAction());
		item.setText("Paste");
		item.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V,
				MenuBarFactory.MENU_MASK));
		menu.add(item);		
		
	}// END: populateMenu
	
	@Override
	public String getMenuName() {
		return "Edit";
	}

	@Override
	public int getPreferredAlignment() {
		return LEFT;
	}

}//END: class
