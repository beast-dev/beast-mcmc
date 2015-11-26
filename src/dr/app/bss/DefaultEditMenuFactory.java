/*
 * DefaultEditMenuFactory.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.bss;

import java.awt.event.KeyEvent;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.text.DefaultEditorKit;

import jam.framework.AbstractFrame;
import jam.framework.MenuBarFactory;
import jam.framework.MenuFactory;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
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
