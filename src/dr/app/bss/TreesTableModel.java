/*
 * TreesTableModel.java
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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.SwingWorker;
import javax.swing.table.AbstractTableModel;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class TreesTableModel extends AbstractTableModel {

	private PartitionDataList dataList;
	private MainFrame frame;

	public final static int TREE_FILE_INDEX = 0;
	public final static int TAXA_SET_INDEX = 1;
	public final static int TAXA_COUNT_INDEX = 2;

	public static String[] COLUMN_NAMES = { "Tree File", "Taxa Set",
			"Taxa count" };

	private static final Class<?>[] COLUMN_TYPES = new Class<?>[] {
			JButton.class, // Tree File
			JButton.class, // Taxa Set
			Integer.class // Taxa
	};

	private TaxaEditor taxaEditor;

	public TreesTableModel(PartitionDataList dataList, MainFrame frame) {
		this.dataList = dataList;
		this.frame = frame;

		addDefaultRow();

	}// END: Constructor

	public void addDefaultRow() {
		dataList.recordsList.add(new TreesTableRecord());
		fireTableDataChanged();
	}// END: addDefaultRow

	public void deleteRow(int row) {

		// remove taxa connected to this row from Taxa panel
		String value = dataList.recordsList.get(row).getName();
		
		Utils.removeTaxaWithAttributeValue(dataList, Utils.TREE_FILENAME, value);

		dataList.recordsList.remove(row);
		fireTableDataChanged();
	}// END: deleteRow

	public void setRow(int row, TreesTableRecord record) {

		// remove taxa connected to this row from Taxa panel
		String value = dataList.recordsList.get(row).getName();
		Utils.removeTaxaWithAttributeValue(dataList, Utils.TREE_FILENAME, value);

		// add taxa with attributes
		if(record.isTreeSet()) {
		applyTaxa(record.getTree());
		}
		
		dataList.recordsList.set(row, record);
		fireTableDataChanged();

	}// END: setRow

	public void applyTaxa(Tree tree) {

		// set attributes to be parsed later in Taxa panel
		Taxa taxa = new Taxa();
		for (Taxon taxon : tree.asList()) {

			double absoluteHeight = Utils.getAbsoluteTaxonHeight(taxon, tree);
			taxon.setAttribute(Utils.ABSOLUTE_HEIGHT, absoluteHeight);
			taxon.setAttribute(Utils.TREE_FILENAME, tree.getId());
			taxa.addTaxon(taxon);
			dataList.allTaxa.addTaxon(taxon);

		}// END: taxon loop

	}// END: applyTaxa

	@Override
	public int getColumnCount() {
		return COLUMN_NAMES.length;
	}

	@Override
	public int getRowCount() {
		return dataList.recordsList.size();
	}

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		return COLUMN_TYPES[columnIndex];
	}// END: getColumnClass

	public boolean isCellEditable(int row, int column) {
		switch (column) {
		case TREE_FILE_INDEX:
			return false;
		case TAXA_COUNT_INDEX:
			return false;
		case TAXA_SET_INDEX:
			return false;
		default:
			return false;
		}
	}// END: isCellEditable

	public String getColumnName(int column) {
		return COLUMN_NAMES[column];
	}// END: getColumnName

	@Override
	public Object getValueAt(int row, int column) {
		switch (column) {

		case TREE_FILE_INDEX:
			JButton treeFileButton = new JButton(Utils.CHOOSE_FILE);
			treeFileButton.addActionListener(new ListenLoadTreeFile(row));
			return treeFileButton;

		case TAXA_SET_INDEX:
			
			JButton taxaEditorButton = new JButton(Utils.EDIT_TAXA_SET);
			taxaEditorButton.addActionListener(new ListenOpenTaxaEditor(row));
			return taxaEditorButton;

		case TAXA_COUNT_INDEX:
			return dataList.recordsList.get(row).getTaxaCount();

		default:
			return "Error";
		}
	}// END: getValueAt

	private class ListenOpenTaxaEditor implements ActionListener {
		private int row;

		public ListenOpenTaxaEditor(int row) {
			this.row = row;
		}

		public void actionPerformed(ActionEvent ev) {

			try {

				taxaEditor = new TaxaEditor(frame, dataList, row);
				taxaEditor.launch();

			} catch (Exception e) {
				Utils.handleException(e);
			}

		}// END: actionPerformed
	}// END: ListenOpenSiteRateModelEditor

	private class ListenLoadTreeFile implements ActionListener {

		private int row;

		public ListenLoadTreeFile(int row) {
			this.row = row;
		}// END: Constructor

		public void actionPerformed(ActionEvent ev) {

			doLoadTreeFile(row);

		}// END: actionPerformed
	}// END: ListenLoadTreeFile

	private void doLoadTreeFile(int row) {

		try {

			JFileChooser chooser = new JFileChooser();
			chooser.setDialogTitle("Select trees file...");
			chooser.setMultiSelectionEnabled(false);
			chooser.setCurrentDirectory(frame.getWorkingDirectory());

			int returnValue = chooser.showOpenDialog(Utils.getActiveFrame());

			if (returnValue == JFileChooser.APPROVE_OPTION) {

				File file = chooser.getSelectedFile();

				if (file != null) {

					loadTreeFile(file, row);

					File tmpDir = chooser.getCurrentDirectory();
					if (tmpDir != null) {
						frame.setWorkingDirectory(tmpDir);
					}

				}// END: file opened check
			}// END: dialog cancelled check

		} catch (Exception e) {
			Utils.handleException(e);
		}// END: try-catch block

	}// END: doLoadTreeFile


	private void loadTreeFile(final File file, final int row) {

		frame.setBusy();
		SwingWorker<Void, Void> worker = new SwingWorker<Void, Void>() {

			// Executed in background thread
			public Void doInBackground() {

				try {

					Tree tree = Utils.importTreeFromFile(file);
					tree.setId(file.getName());

					setRow(row, new TreesTableRecord(tree.getId(), tree));

				} catch (Exception e) {
					Utils.handleException(e);
				}// END: try-catch block

				return null;
			}// END: doInBackground()

			// Executed in event dispatch thread
			public void done() {
				frame.setIdle();
				frame.fireTaxaChanged();
			}// END: done
		};

		worker.execute();

	}// END: loadTreeFile

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
