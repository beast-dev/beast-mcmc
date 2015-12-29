/*
 * GeneralPreferencesSection.java
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

package dr.app.mapper.application;

import jam.preferences.PreferencesSection;
import jam.util.IconUtils;

import javax.swing.*;
import java.util.prefs.Preferences;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class GeneralPreferencesSection implements PreferencesSection {
	Icon projectToolIcon = IconUtils.getIcon(MapperApp.class, "images/prefsGeneral.png");

	public String getTitle() {
		return "General";
	}

	public Icon getIcon() {
		return projectToolIcon;
	}

	public JPanel getPanel() {
		JPanel panel = new JPanel();
		panel.add(generalCheck);
		return panel;
	}

	public void retrievePreferences() {
		Preferences prefs = Preferences.userNodeForPackage(MapperApp.class);
		generalCheck.setSelected(prefs.getBoolean("general_check", true));
	}

	public void storePreferences() {
		Preferences prefs = Preferences.userNodeForPackage(MapperApp.class);
		prefs.putBoolean("general_check", generalCheck.isSelected());
	}

	JCheckBox generalCheck = new JCheckBox("The preferences window is not implemented yet.");
}
