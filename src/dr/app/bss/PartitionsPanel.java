/*
 * PartitionsPanel.java
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
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableColumn;


/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class PartitionsPanel extends JPanel implements Exportable {

	private PartitionDataList dataList = null;

	private JTable partitionTable = null;
	private PartitionTableModel partitionTableModel = null;
	private JTableComboBoxCellEditor jTableComboBoxCellEditor = null;
	
	private TableColumnHider hider;
	private JScrollPane scrollPane;
	private TableColumn column;
	private int columnIndex;
	private int partitionsCount;

	private Action addPartitionAction = new AbstractAction("+") {
		public void actionPerformed(ActionEvent ae) {

			// partitionTableModel.copyPreviousRow();
			partitionTableModel.addDefaultRow();
			setPartitions();
			
		}// END: actionPerformed
	};

	private Action removePartitionAction = new AbstractAction("-") {
		public void actionPerformed(ActionEvent ae) {
			if (partitionsCount > 1) {
				
				partitionTableModel.deleteRow(partitionsCount - 1);
				setPartitions();
			
			}
		}// END: actionPerformed
	};

	public PartitionsPanel(PartitionDataList dataList) {

		this.dataList = dataList;

		partitionTable = new JTable();
		partitionTable.getTableHeader().setReorderingAllowed(false);
		partitionTable.addMouseListener(new JTableButtonMouseListener(
				partitionTable));

		partitionTableModel = new PartitionTableModel(this.dataList);
		partitionTableModel
				.addTableModelListener(new PartitionTableModelListener());
		partitionTable.setModel(partitionTableModel);

		hider = new TableColumnHider(partitionTable);

		setLayout(new BorderLayout());

		scrollPane = new JScrollPane(partitionTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		RowNumberTable rowNumberTable = new RowNumberTable(partitionTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER,
				rowNumberTable.getTableHeader());

		scrollPane.getViewport().setOpaque(false);

		add(scrollPane, BorderLayout.CENTER);

		//TODO: Input
		columnIndex = PartitionTableModel.DATA_INDEX;
		column = partitionTable.getColumnModel().getColumn(
				columnIndex);
		column.setCellRenderer(new JTableComboBoxCellRenderer(columnIndex));
	    jTableComboBoxCellEditor = new JTableComboBoxCellEditor(this.dataList); 
		column.setCellEditor(jTableComboBoxCellEditor);
		column.setMinWidth(100);
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.FROM_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		column.setPreferredWidth(80);
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.TO_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		column.setPreferredWidth(80);
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.EVERY_INDEX);
		column.setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 2,
						0, 2)));
		column.setPreferredWidth(80);

		columnIndex = PartitionTableModel.DATA_TYPE_INDEX;
		column = partitionTable.getColumnModel().getColumn(
				columnIndex);
		column.setCellEditor(new JTableComboBoxCellEditor(this.dataList));
		column.setCellRenderer(new JTableComboBoxCellRenderer(columnIndex));
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.DEMOGRAPHIC_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.BRANCH_SUBSTITUTION_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.FREQUENCY_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());
		
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.SITE_RATE_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.CLOCK_RATE_MODEL_INDEX);
		column.setCellRenderer(new JTableButtonCellRenderer());
		column.setCellEditor(new JTableButtonCellEditor());

		String toolTipText = "Leave this field empty to start from a random root sequence.";
		column = partitionTable.getColumnModel().getColumn(
				PartitionTableModel.ROOT_SEQUENCE_INDEX);
		JTableButtonCellRenderer jTableButtonCellRenderer = new JTableButtonCellRenderer();
		jTableButtonCellRenderer.setColumnToolTipText(toolTipText);
		column.setCellRenderer(jTableButtonCellRenderer);
		column.setCellEditor(new JTableButtonCellEditor());
		
		ActionPanel actionPanel = new ActionPanel(false);
		actionPanel.setAddAction(addPartitionAction);
		actionPanel.setRemoveAction(removePartitionAction);
		add(actionPanel, BorderLayout.SOUTH);

		setPartitions();

	}// END: Constructor

	private void setPartitions() {

		partitionsCount = dataList.size();

		addPartitionAction.setEnabled(true);
		if (partitionsCount == 1) {
			removePartitionAction.setEnabled(false);
		} else {
			removePartitionAction.setEnabled(true);
		}

		ColumnResizer.adjustColumnPreferredWidths(partitionTable);
	}// END: setPartitions

	// Listen to tree choices, set tree model in partition data
	private class PartitionTableModelListener implements TableModelListener {

		public void tableChanged(TableModelEvent ev) {

			if (ev.getType() == TableModelEvent.UPDATE) {
				int row = ev.getFirstRow();
				int column = ev.getColumn();

				if (column == PartitionTableModel.DATA_INDEX) {

					TreesTableRecord value = (TreesTableRecord) partitionTableModel.getValueAt(row,
							column);
					
					dataList.get(row).record = (TreesTableRecord) value;

				}// END: column check

			}// END: event check

//			 frame.collectAllSettings();

		}// END: tableChanged
	}// END: InteractiveTableModelListener

	public void hideTreeColumn() {
		hider.hide(PartitionTableModel.COLUMN_NAMES[PartitionTableModel.DATA_INDEX]);
	}

	public void showTreeColumn() {
		hider.show(PartitionTableModel.COLUMN_NAMES[PartitionTableModel.DATA_INDEX]);
	}

	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	public void updatePartitionTable(PartitionDataList dataList) {
		partitionTableModel.setDataList(dataList);
		setDataList(dataList);
		setPartitions();
		jTableComboBoxCellEditor.setDataList(dataList);
		partitionTableModel.fireTableDataChanged();
	}// END: updatePartitionTable

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}

}// END: class
