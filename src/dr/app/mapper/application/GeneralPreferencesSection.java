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
