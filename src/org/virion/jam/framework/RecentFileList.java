/*
 * RecentFileList.java
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
import javax.swing.event.EventListenerList;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Properties;
import java.util.Vector;

/**
 * A class for maintaining a "Recent File List".
 * The recent file list can be stored between program invocations in
 * a properties file. One or more RecentFileLists can easily
 * be embedded in a JMenu.
 *
 * @author Tony Johnson (tonyj@slac.stanford.edu)
 * @version $Id: RecentFileList.java,v 1.2 2006/09/09 18:16:16 rambaut Exp $
 */
public class RecentFileList implements ActionListener {
    /**
     * Create a RecentFileList
     */
    public RecentFileList() {
        this(null);
    }

    /**
     * Create a RecentFileList with a given maximum length
     *
     * @param size the maximum number of files to remember
     */
    public RecentFileList(int size) {
        this(null, size);
    }

    /**
     * Create a recent file list. The type parameter is used to
     * prefix entries in the properties file, so that multiple
     * RecentFileLists can be used in an application.
     *
     * @param type The prefix to use
     */
    public RecentFileList(String type) {
        this(type, 4);
    }

    /**
     * Create a recent file list with a given type and size
     *
     * @param type The prefix to use
     * @param size the maximum number of files to remember
     */
    public RecentFileList(String type, int size) {
        files = new Vector(size);
        this.size = size;
        this.type = type;
        used = 0;
    }

    public void addActionListener(ActionListener l) {
        listenerList.add(ActionListener.class, l);
    }

    public void removeActionListener(ActionListener l) {
        listenerList.remove(ActionListener.class, l);
    }

    /**
     * An action event is fired when the user selects one of the
     * files from a menu. The "actionCommand" in the event will be
     * set to the name of the selected file.
     */
    protected void fireActionPerformed(ActionEvent e) {
        // Guaranteed to return a non-null array
        Object[] listeners = listenerList.getListenerList();
        // Process the listeners last to first, notifying
        // those that are interested in this event
        for (int i = listeners.length - 2; i >= 0; i -= 2) {
            if (listeners[i] == ActionListener.class) {
                ((ActionListener) listeners[i + 1]).actionPerformed(e);
            }
        }
    }

    public void actionPerformed(ActionEvent e) {
        fireActionPerformed(e);
    }

    /**
     * Save the recent file list in a Properties set
     *
     * @param props The Properties set to save the files in
     */
    public void save(Properties props) {
        String key = "RecentFile_" + (type != null ? type + "_" : "");
        for (int i = 0; i < used; i++) {
            props.put(key + i, files.elementAt(i).toString());
        }
        props.put(key + "Used", String.valueOf(used));
    }

    /**
     * Load the recent file list from a Properties set
     *
     * @param props The Properties set to load from
     */
    public void load(Properties props) {
        String key = "RecentFile_" + (type != null ? type + "_" : "");
        files.removeAllElements();
        used = Integer.parseInt(props.getProperty(key + "Used", "0"));
        for (int i = 0; i < used; i++) {
            String value = props.getProperty(key + i);
            if (value == null) break;
            files.addElement(value);
        }
    }

    /**
     * Add a file to the list
     *
     * @param f The file to add
     */
    public void add(File f) {
        try {
            add(f.getCanonicalPath());
        } catch (java.io.IOException x) {
            add(f.getAbsolutePath());
        }
    }

    /**
     * Remove a file from the list
     *
     * @param f Remove a file from the list
     */
    public void remove(File f) {
        remove(f.getName());
    }

    /**
     * Add a file to the list
     *
     * @param name The name of the file to add
     */
    public void add(String name) {
        int pos = files.indexOf(name);
        if (pos > 0) {
            files.removeElementAt(pos);
            files.insertElementAt(name, 0);
        } else if (pos != 0) {
            if (used == size)
                files.removeElementAt(size - 1);
            else
                used++;
            files.insertElementAt(name, 0);
        }
    }

    /**
     * Remove a file from the list
     *
     * @param name The name of the file to remove
     */
    public void remove(String name) {
        if (files.removeElement(name)) used--;
    }

    /**
     * Adds the recent file list to a menu.
     * The files will be added at the end of the menu, with a
     * separator before the files (if there are >0 files in the list)
     */
    public void buildMenu(JMenu menu) {
        if (used > 0) {
            menu.addSeparator();
            for (int i = 0; i < used; i++) {
                JMenuItem item = new JMenuItem(String.valueOf(i + 1) + " " + files.elementAt(i));
                item.setActionCommand(files.elementAt(i).toString());
                if (size < 9) item.setMnemonic(Character.forDigit(i + 1, 10));
                item.addActionListener(this);
                menu.add(item);
            }
        }
    }

    private Vector files;
    private String type;
    private int size;
    private int used;
    private EventListenerList listenerList = new EventListenerList();
}
