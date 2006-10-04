/*
 * DefaultMenuBarFactory.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.framework;

import javax.swing.*;
import java.util.*;

/**
 * @author rambaut
 *         Date: Dec 26, 2004
 *         Time: 10:55:55 AM
 */
public class DefaultMenuBarFactory implements MenuBarFactory {

	private final List menuFactories = new ArrayList();
    private final List permanentMenuFactories = new ArrayList();
    private boolean populatedMenu = false;
    private JMenuBar menuBar = null;
    AbstractFrame frame = null;

	public final void populateMenuBar(JMenuBar menuBar, AbstractFrame frame) {
        this.menuBar = menuBar;
        this.frame = frame;
		Map menus = new HashMap();
		Map order = new TreeMap();
        List l = new ArrayList();
        l.addAll(menuFactories);
        l.addAll(permanentMenuFactories);
        Iterator iter = l.iterator();

		int leftOrder = 0;
		int centerOrder = 1000;
		int rightOrder = 10000;

		while (iter.hasNext()) {
			MenuFactory menuFactory = (MenuFactory)iter.next();
			String name = menuFactory.getMenuName();
			JMenu menu = (JMenu)menus.get(name);
			if (menu == null) {
				int alignment = menuFactory.getPreferredAlignment();
				menu = new JMenu(name);
				menus.put(name, menu);
				switch (alignment) {
					case MenuFactory.LEFT:
						order.put(new Integer(leftOrder), name);
						leftOrder++;
					break;
					case MenuFactory.CENTER:
						order.put(new Integer(centerOrder), name);
						centerOrder++;
					break;
					case MenuFactory.RIGHT:
						order.put(new Integer(rightOrder), name);
						rightOrder--;
					break;
				}
			}

			menuFactory.populateMenu(menu, frame);
		}

//        Iterator iter3 = menuFactories.iterator();
//
//        while (iter3.hasNext()) {
//            MenuFactory menuFactory = (MenuFactory)iter3.next();
//            String name = menuFactory.getMenuName();
//            JMenu menu = (JMenu)menus.get(name);
//            if (menu == null) {
//                int alignment = menuFactory.getPreferredAlignment();
//                menu = new JMenu(name);
//                menus.put(name, menu);
//                switch (alignment) {
//                    case MenuFactory.LEFT:
//                        order.put(new Integer(leftOrder), name);
//                        leftOrder++;
//                    break;
//                    case MenuFactory.CENTER:
//                        order.put(new Integer(centerOrder), name);
//                        centerOrder++;
//                    break;
//                    case MenuFactory.RIGHT:
//                        order.put(new Integer(rightOrder), name);
//                        rightOrder--;
//                    break;
//                }
//            }
//
//            menuFactory.populateMenu(menu, frame);
//        }

		Iterator iter2 = order.keySet().iterator();
		while (iter2.hasNext()) {
			Integer i = (Integer)iter2.next();
			String name = (String)order.get(i);
			JMenu menu = (JMenu)menus.get(name);

			menuBar.add(menu);
		}
        populatedMenu = true;
    }

    public final void deregisterMenuFactories() {
        menuFactories.removeAll(menuFactories);
	}

    public final void registerPermanentMenuFactory(MenuFactory menuFactory) {
        permanentMenuFactories.add(menuFactory);
        if(populatedMenu) {
            menuBar.removeAll();
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    populateMenuBar(menuBar, frame);
                    frame.setJMenuBar(menuBar);
                }
            });
        }
    }

	public final void registerMenuFactory(MenuFactory menuFactory) {
		menuFactories.add(menuFactory);
        if(populatedMenu) {
            menuBar.removeAll();
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    populateMenuBar(menuBar, frame);
                    frame.setJMenuBar(menuBar);
	}
            });
        }
    }

}
