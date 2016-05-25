/*
 * TaxonListTable.java
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

package dr.app.gui.table;

import dr.evolution.util.*;
import jam.table.*;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.Iterator;
import java.util.TreeSet;

public class TaxonListTable extends JTable implements MutableTaxonListListener {

	/**
	 *
	 */
	private static final long serialVersionUID = -5706027520036473157L;

	final static String TAXON_NAME = "Taxon";

	TaxaTableModel taxaTableModel = null;
	TaxonList taxonList = null;

	public TaxonListTable(TaxonList taxonList) {

		this.taxonList = taxonList;

		if (taxonList instanceof MutableTaxonList) {
			((MutableTaxonList)taxonList).addMutableTaxonListListener(this);
		}

		String[] names = getColumnNames();

		// Create a model of the data.
		taxaTableModel = new TaxaTableModel(names, taxonList);

		// Create the table
		//JTable tableView = new JTable(dataModel);
		TableSorter sorter = new TableSorter(taxaTableModel); //ADDED THIS
		setModel(sorter);
		sorter.addTableModelListener(this); //ADDED THIS

		// Turn off auto-resizing so that we can set column sizes programmatically.
		// In this mode, all columns will get their preferred widths, as set below.
		setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		setColumnSelectionAllowed(true);

		// Set the widths of the columns.
		getColumn(TAXON_NAME).setPreferredWidth(120);
		for (int i = 0; i < getColumnCount(); i++) {
			getColumnModel().getColumn(i).setCellRenderer(
				new TableRenderer(true, SwingConstants.LEFT, new Insets(0, 4, 0, 0)));
		}

		setDefaultEditor(dr.evolution.util.Date.class, new DateCellEditor());
	}

	public void updateColumns() {
		taxaTableModel.setColumnNames(getColumnNames());
	}

	private String[] getColumnNames() {
		TreeSet attributes = new TreeSet();

		for (int i = 0; i < taxonList.getTaxonCount(); i++) {
			Taxon taxon = taxonList.getTaxon(i);
			Iterator iter = taxon.getAttributeNames();
			while (iter.hasNext()) {
				Object attribute = iter.next();
				attributes.add(attribute);
			}

		}

		String[] names = new String[attributes.size() + 1];
		names[0] = "Taxon";
		Iterator iter = attributes.iterator();
		int i = 1;
		while (iter.hasNext()) {
			names[i] = iter.next().toString();
			i++;
		}

		return names;
	}

	//********************************************************************
	// MutableTaxonListListener interface
	//********************************************************************

	public void taxonAdded(TaxonList taxonList, Taxon taxon) {

		taxaTableModel.fireTableDataChanged();
	}

	public void taxonRemoved(TaxonList taxonList, Taxon taxon) {
		taxaTableModel.fireTableDataChanged();
	}

	public void taxaChanged(TaxonList taxonList) {
		updateColumns();
	}

	//************************************************************************
	// TableModel inner class
	//************************************************************************

	class TaxaTableModel extends AbstractTableModel {

		/**
		 *
		 */
		private static final long serialVersionUID = 367162013042095775L;
		String[] names = null;
		TaxonList taxonList;

		public TaxaTableModel(String[] names, TaxonList taxonList) {
			this.names = names;
			this.taxonList = taxonList;
		}

		public void setColumnNames(String[] names) {
			this.names = names;
			fireTableStructureChanged();
		}

		// These methods always need to be implemented.
		public int getColumnCount() { return names.length; }
		public int getRowCount() { return taxonList.getTaxonCount();}

		public Object getValueAt(int row, int column) {

			if (column == 0) {
				return taxonList.getTaxon(row).getId();
			} else {
				return taxonList.getTaxonAttribute(row, names[column]);
			}
		}

		// The default implementations of these methods in
		// AbstractTableModel would work, but we can refine them.
		public String getColumnName(int column) {return names[column];}
		public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
		public boolean isCellEditable(int row, int column) {

			if (taxonList instanceof MutableTaxonList) {
				if (column == 0) return true;

				Object item = taxonList.getTaxonAttribute(row, names[column]);
				if (item instanceof String) return true;
				if (item instanceof dr.evolution.util.Date) return true;
			}
			return false;
		}

		public void setValueAt(Object aValue, int row, int column) {

			if (taxonList instanceof MutableTaxonList) {

				if (column == 0) {
					((MutableTaxonList)taxonList).setTaxonId(row, (String)aValue);
				} else {
					((MutableTaxonList)taxonList).setTaxonAttribute(row, names[column], aValue);
				}

			}
		}

	}
}
