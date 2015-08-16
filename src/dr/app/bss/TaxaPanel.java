/*
 * TaxaPanel.java
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
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import java.awt.BorderLayout;
import java.awt.Insets;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.TableColumn;

import dr.app.gui.table.TableEditorStopper;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class TaxaPanel extends JPanel implements Exportable {

	private PartitionDataList dataList = null;

	private JScrollPane scrollPane = new JScrollPane();
	private JTable taxaTable = null;
	private TaxaTableModel taxaTableModel = null;
	private TableColumn column;

	public TaxaPanel(PartitionDataList dataList) {

		this.dataList = dataList;

		taxaTable = new JTable();

		taxaTableModel = new TaxaTableModel(this.dataList);
		taxaTable.setModel(taxaTableModel);

		setLayout(new BorderLayout());

		taxaTable.getTableHeader().setReorderingAllowed(false);

		taxaTable.getTableHeader()
				.setDefaultRenderer(
						new HeaderRenderer(SwingConstants.LEFT, new Insets(0,
								2, 0, 2)));

		column = taxaTable.getColumnModel()
				.getColumn(TaxaTableModel.NAME_INDEX);
		column.setCellRenderer(new TableRenderer(SwingConstants.LEFT,
				new Insets(0, 2, 0, 2)));
		column.setPreferredWidth(80);

		column = taxaTable.getColumnModel().getColumn(
				TaxaTableModel.HEIGHT_INDEX);
		column.setCellRenderer(new TableRenderer(SwingConstants.LEFT,
				new Insets(0, 2, 0, 2)));
		column.setPreferredWidth(80);

		column = taxaTable.getColumnModel().getColumn(
				TaxaTableModel.TAXA_SET_INDEX);
		column.setCellRenderer(new TableRenderer(SwingConstants.LEFT,
				new Insets(0, 2, 0, 2)));
		column.setPreferredWidth(80);

		TableEditorStopper.ensureEditingStopWhenTableLosesFocus(taxaTable);

		scrollPane = new JScrollPane(taxaTable,
				JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		RowNumberTable rowNumberTable = new RowNumberTable(taxaTable);
		scrollPane.setRowHeaderView(rowNumberTable);
		scrollPane.setCorner(JScrollPane.UPPER_LEFT_CORNER,
				rowNumberTable.getTableHeader());

		scrollPane.getViewport().setOpaque(false);

		setOpaque(false);
		setBorder(new BorderUIResource.EmptyBorderUIResource(
				new java.awt.Insets(12, 12, 12, 12)));
		setLayout(new BorderLayout(0, 0));

		add(scrollPane, BorderLayout.CENTER);

	}// END: Constructor

	public JComponent getExportableComponent() {
		return taxaTable;
	}

	public void updateTaxaTable(PartitionDataList dataList) {
		taxaTableModel.setDataList(dataList);
		setDataList(dataList);
		fireTaxaChanged();
	}// END: updateTaxaTable

	public void fireTaxaChanged() {
		taxaTableModel.fireTaxaChanged();
	}// END: fireTableDataChanged

	public void setDataList(PartitionDataList dataList) {
		this.dataList = dataList;
	}// END: setDataList

}// END: class
