/*
 * TreeStatPanel.java
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

package dr.app.treestat;

import jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that displays information about trees
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TreeStatPanel.java,v 1.3 2006/05/09 10:24:27 rambaut Exp $
 */
public class TreeStatPanel extends javax.swing.JPanel implements Exportable {

	/**
	 *
	 */
	private static final long serialVersionUID = 2437334458007083790L;
	private JTabbedPane tabbedPane = new JTabbedPane();
	private TaxonSetsPanel taxonSetsPanel;
	//private CharactersPanel charactersPanel;
	private StatisticsPanel statisticsPanel;

//	final Icon traceIcon = new ImageIcon(Utils.getImage(this, "/images/trace-small-icon.gif"));

	/** Creates new form TreeStatPanel */
	public TreeStatPanel(TreeStatFrame frame, TreeStatData treeStatData) {
		taxonSetsPanel = new TaxonSetsPanel(frame, treeStatData);
		//charactersPanel = new CharactersPanel(frame, treeStatData);
		statisticsPanel = new StatisticsPanel(frame, treeStatData);

		tabbedPane.addTab("Taxon Sets", null, taxonSetsPanel);
		//tabbedPane.addTab("Characters", null, charactersPanel);
		tabbedPane.addTab("Statistics", null, statisticsPanel);

		setLayout(new BorderLayout());
		add(tabbedPane, BorderLayout.CENTER);
	}

	/** This function calls dataChanged for each panel */
	public void fireDataChanged() {

		taxonSetsPanel.dataChanged();
		//charactersPanel.dataChanged();
		statisticsPanel.dataChanged();
	}

    public void doCopy() {
    	//summaryPanel.copyToClipboard();
    	/*switch (tabbedPane.getSelectedIndex()) {
			case 0: summaryPanel.copyToClipboard(); break;
			case 1: densityPanel.copyToClipboard(); break;
			case 2: tracePanel.copyToClipboard(); break;
		}*/
	}

    public JComponent getExportableComponent() {

		JComponent exportable = null;
		Component comp = tabbedPane.getSelectedComponent();

		if (comp instanceof Exportable) {
			exportable = ((Exportable)comp).getExportableComponent();
		} else if (comp instanceof JComponent) {
			exportable = (JComponent)comp;
		}

		return exportable;
	}

	//************************************************************************
	// private methods
	//************************************************************************

}
