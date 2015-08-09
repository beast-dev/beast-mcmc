/*
 * TreesPanel.java
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

import jam.framework.Exportable;
import jam.panels.ActionPanel;
import jam.table.TableRenderer;

import java.awt.BorderLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.TableColumn;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class TreesPanel extends JPanel implements Exportable {

	private PartitionDataList dataList;
	private MainFrame frame;

	private JTable treesTable = null;
	private TreesTableModel treesTableModel = null;

	private JScrollPane scrollPane;
	private TableColumn column;
	private int rowCount;

	private Action addTreeAction = new AbstractAction("+") {
		public void actionPerformed(ActionEvent ae) {

			treesTableModel.addDefaultRow();
			setTrees();

		}// END: actionPerformed
	};

	private Action removeTreeAction = new AbstractAction("-") {
		public void actionPerformed(ActionEvent ae) {
			if (rowCount > 1) {

				treesTableModel.deleteRow(rowCount - 1);
				frame.fireTaxaChanged();
				setTrees();

			}
		}// END: actionPerformed
	};

	public TreesPanel(MainFrame frame, PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		treesTable = new JTable();
		treesTable.getTableHeader().setReorderingAllowed(false);
		treesTable.addMouseListener(new JTableButtonMouseListener(treesTable));

		treesTableModel = new TreesTableModel(this.dataList, this.frame);
		treesTable.setModel(treesTableModel);

		setLayout(new BorderLayout());

		scrollPane = new JScrollPane(treesTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		RowNumberTable rowNumberTable = new RowNumberTable(treesTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER,
				rowNumberTable.getTableHeader());
		scrollPane.getViewport().setOpaque(false);
		add(scrollPane, BorderLayout.CENTER);

		setTreesColumn(this.dataList);
		setTaxaSetColumn(this.dataList);

		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TAXA_COUNT_INDEX);
		column.setCellRenderer(new TableRenderer(SwingConstants.LEFT,
				new Insets(0, 2, 0, 2)));

		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addTreeAction);
		actionPanel.setRemoveAction(removeTreeAction);
		add(actionPanel, BorderLayout.SOUTH);

		setTrees();

	}// END: Constructor

	private void setTaxaSetColumn(PartitionDataList dataList) {
		
		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TAXA_SET_INDEX);
		// pass dataList for labels on buttons
		column.setCellRenderer(new JTableButtonCellRenderer(
				dataList,TreesTableModel.TAXA_SET_INDEX
				));
		column.setCellEditor(new JTableButtonCellEditor());
		
	}
	
	private void setTreesColumn(PartitionDataList dataList) {

		column = treesTable.getColumnModel().getColumn(
				TreesTableModel.TREE_FILE_INDEX);
		// pass dataList for labels on buttons
		column.setCellRenderer(new JTableButtonCellRenderer(dataList, TreesTableModel.TREE_FILE_INDEX));
		column.setCellEditor(new JTableButtonCellEditor());

	}// END: setTreesColumn

	private void setTrees() {

		rowCount = dataList.recordsList.size();

		addTreeAction.setEnabled(true);
		if (rowCount == 1) {
			removeTreeAction.setEnabled(false);
		} else {
			removeTreeAction.setEnabled(true);
		}

		ColumnResizer.adjustColumnPreferredWidths(treesTable);
	}// END: setPartitions

	public void updateTreesTable(PartitionDataList dataList) {
		treesTableModel.setDataList(dataList);
		setDataList(dataList);
		setTreesColumn(dataList);
		setTaxaSetColumn(dataList);
		setTrees();
		fireTableDataChanged();
	}// END: updateTreesTable

	public void fireTableDataChanged() {
		treesTableModel.fireTableDataChanged();
	}
	
	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
