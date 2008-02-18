/*
 * DataPanel.java
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

package dr.app.coalgen;

import dr.evolution.util.Date;
import dr.evolution.util.Units;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.DateCellEditor;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableRenderer;
import org.virion.jam.table.TableSorter;
import org.virion.jam.components.RealNumberField;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id$
 */
public class DataPanel extends JPanel implements Exportable {

	JScrollPane scrollPane = new JScrollPane();
	JTable dataTable = null;
	DataTableModel dataTableModel = null;

	ClearDatesAction clearDatesAction = new ClearDatesAction();
	GuessDatesAction guessDatesAction = new GuessDatesAction();
	
	JComboBox unitsCombo = new JComboBox(new String[] {"Years", "Months", "Days"});
	JComboBox directionCombo = new JComboBox(new String[] {"Before the present", "Since some time in the past"});
	//RealNumberField originField = new RealNumberField(0.0, Double.POSITIVE_INFINITY);
	
	JComboBox translationCombo = new JComboBox();

	TableRenderer sequenceRenderer = null;
		
	CoalGenFrame frame = null;
	
	CoalGenData data = null;
	
	public DataPanel(CoalGenFrame frame, CoalGenData data) {
	
		this.frame = frame;
		this.data = data;
		
		dataTableModel = new DataTableModel();
        TableSorter sorter = new TableSorter(dataTableModel);
		dataTable = new JTable(sorter);
	
        sorter.setTableHeader(dataTable.getTableHeader());

		dataTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
		dataTable.getTableHeader().setReorderingAllowed(false); 
		dataTable.getTableHeader().setDefaultRenderer(
			new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4))); 

		dataTable.getColumnModel().getColumn(0).setCellRenderer(
			new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
		dataTable.getColumnModel().getColumn(0).setPreferredWidth(80);
		dataTable.getColumnModel().getColumn(1).setCellRenderer(
			new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4))); 
		dataTable.getColumnModel().getColumn(1).setPreferredWidth(80);
		dataTable.getColumnModel().getColumn(1).setCellEditor(
			new DateCellEditor());

		dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent evt) { selectionChanged(); }
		});

   		scrollPane = new JScrollPane(dataTable, 
										JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, 
										JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		scrollPane.setOpaque(false);

		clearDatesAction.setEnabled(false);
		guessDatesAction.setEnabled(false);
		unitsCombo.setOpaque(false);
		unitsCombo.setEnabled(false);
		unitsCombo.setFont(unitsCombo.getFont().deriveFont(12.0f));
		directionCombo.setOpaque(false);
		directionCombo.setEnabled(false);
		directionCombo.setFont(directionCombo.getFont().deriveFont(12.0f));
		//originField.setEnabled(false);
		//originField.setValue(0.0);
		//originField.setColumns(12);
		
		JToolBar toolBar1 = new JToolBar();
		toolBar1.setFloatable(false);
		toolBar1.setOpaque(false);
//		toolBar1.setLayout(new BoxLayout(toolBar1, javax.swing.BoxLayout.X_AXIS));
		toolBar1.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
		toolBar1.add(clearDatesAction);
		toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));
		toolBar1.add(guessDatesAction);
		toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));
		toolBar1.add(new JLabel("Dates specified as "));
		toolBar1.add(unitsCombo);
		toolBar1.add(directionCombo);
		//toolBar.add(originField);

		setOpaque(false);
		setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
   		setLayout(new BorderLayout(0,0));
		add(toolBar1, "North");
		add(scrollPane, "Center");

		ItemListener listener =	new ItemListener() {
			public void itemStateChanged(ItemEvent ev) { timeScaleChanged(); } 
		};
		unitsCombo.addItemListener(listener);
		directionCombo.addItemListener(listener);
		//originField.addKeyListener(new java.awt.event.KeyAdapter() {
		//	public void keyTyped(java.awt.event.KeyEvent ev) {
		//		timeScaleChanged();
		//	}});

	}
	
	public final void dataChanged() {
		
		if (data.taxonList.getTaxonCount() > 0) {			
			clearDatesAction.setEnabled(true);
			guessDatesAction.setEnabled(true);
			unitsCombo.setEnabled(true);
			directionCombo.setEnabled(true);

			//originField.setEnabled(true);

		}
		
		dataTableModel.fireTableDataChanged();
	}
	
	public final void timeScaleChanged() {
		int units = Units.YEARS;
		switch (unitsCombo.getSelectedIndex()) {
			case 0: units = Units.YEARS; break;
			case 1: units = Units.MONTHS; break;
			case 2: units = Units.DAYS; break;
		}
		
		boolean backwards = directionCombo.getSelectedIndex() == 0;
		
		//double origin = originField.getValue().doubleValue();
		
		for (int i = 0; i < data.taxonList.getTaxonCount(); i++) {
			Date date = data.taxonList.getTaxon(i).getDate();
			double d = date.getTimeValue();
			
			Date newDate = createDate(d, units, backwards, 0.0);
			
			data.taxonList.getTaxon(i).setDate(newDate);
		}	
						
		dataTableModel.fireTableDataChanged();
		frame.fireDataChanged();
	}
	
	private Date createDate(double timeValue, int units, boolean backwards, double origin) {
		if (backwards) {
			return Date.createTimeAgoFromOrigin(timeValue, units, origin);
		} else {
			return Date.createTimeSinceOrigin(timeValue, units, origin);
		}
	}
		
    public JComponent getExportableComponent() {
		return dataTable;
	} 	
      
	public void selectionChanged() {
		
		int[] selRows = dataTable.getSelectedRows();
		if (selRows == null || selRows.length == 0) {
			frame.dataSelectionChanged(false);
		} else {
			frame.dataSelectionChanged(true);
		}
	}

	public void deleteSelection() {
/*		int option = JOptionPane.showConfirmDialog(this, "Are you sure you wish to delete\n"+
														 "the selected taxa?\n"+
														 "This operation cannot be undone.", 
													"Warning",
													JOptionPane.YES_NO_OPTION,
													JOptionPane.WARNING_MESSAGE);
											
		if (option == JOptionPane.YES_OPTION) {
			int[] selRows = dataTable.getSelectedRows();
			String[] names = new String[selRows.length];
			
			TableModel model = dataTable.getModel();
			
			for (int i = 0; i < names.length; i++) {
				names[i] = (String)model.getValueAt(selRows[i], 0);
			}
			
			for (int i = 0; i < names.length; i++) {
				int index = data.taxonList.getTaxonIndex(names[i]);
				data.taxonList.removeTaxon(index);
			}
			
			if (options.originalAlignment.getTaxonCount() == 0) {
				// if all the sequences are deleted we may as well throw
				// away the alignment...
			
				options.originalAlignment = null;
				options.alignment = null;
			}
			
			dataTableModel.fireTableDataChanged();
			frame.dataChanged();
		}
		*/
	}

  	public void clearDates() {
		for (int i = 0; i < data.taxonList.getTaxonCount(); i++) {
			java.util.Date origin = new java.util.Date(0);
			
			double d = 0.0;
			
			Date date = Date.createTimeSinceOrigin(d, Units.YEARS, origin);
			data.taxonList.getTaxon(i).setAttribute("date", date);
		}
		
		// adjust the dates to the current timescale...	
		timeScaleChanged();

		dataTableModel.fireTableDataChanged();
		frame.fireDataChanged();
	}
	
  	public void guessDates() {
  	
  		OptionsPanel optionPanel = new OptionsPanel();
  		
  		optionPanel.addLabel("The date is given by a numerical field in the taxon label that is:");
  	
		final JLabel orderLabel = new JLabel("Defined by its order:");
		final JComboBox orderCombo = new JComboBox(new String[] {"first", "second", "third",
									"fourth", "fourth from last", 
									"third from last", "second from last", "last"});
								
		optionPanel.addComponents(orderLabel, orderCombo);
		optionPanel.addSeparator();
									
		final JCheckBox prefixCheckBox = new JCheckBox("Defined by a prefix", false);
 		final JTextField prefixText = new JTextField(16);
		prefixText.setEnabled(false);
		optionPanel.addComponents(prefixCheckBox, prefixText);
 		optionPanel.addSeparator();
	
		final JCheckBox offsetCheck = new JCheckBox("Add the following value to each: ", false);
 		final RealNumberField offsetText = new RealNumberField();
 		offsetText.setValue(1900);
 		offsetText.setColumns(16);
		offsetText.setEnabled(false);
		offsetCheck.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				offsetText.setEnabled(offsetCheck.isSelected());
			}
        });
		optionPanel.addComponents(offsetCheck, offsetText);

		final JCheckBox unlessCheck = new JCheckBox("...unless less than:", false);
 		final RealNumberField unlessText = new RealNumberField();
 		unlessText.setValue(4);
 		unlessText.setColumns(16);
		unlessText.setEnabled(false);
		optionPanel.addComponents(unlessCheck, unlessText);

 		final RealNumberField offset2Text = new RealNumberField();
 		offset2Text.setValue(2000);
 		offset2Text.setColumns(16);
		offset2Text.setEnabled(false);
		optionPanel.addComponentWithLabel("...in which case add:", offset2Text);

		unlessCheck.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				unlessText.setEnabled(unlessCheck.isSelected());
				offset2Text.setEnabled(unlessCheck.isSelected());
			}
        });

        prefixCheckBox.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				prefixText.setEnabled(prefixCheckBox.isSelected());
			}
        });
  		
		JOptionPane optionPane = new JOptionPane(optionPanel,
                                    JOptionPane.QUESTION_MESSAGE,
                                    JOptionPane.OK_CANCEL_OPTION,
                                    null,
                                    null,
                                    null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		JDialog dialog = optionPane.createDialog(frame, "Guess Dates");	
//		dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
   		dialog.pack();
		dialog.show();
		
		if (optionPane.getValue() == null) {
			return;
		}
		
		int value = ((Integer)optionPane.getValue()).intValue();
		if (value == -1 || value == JOptionPane.CANCEL_OPTION) {
		    return;
		}
		
		for (int i = 0; i < data.taxonList.getTaxonCount(); i++) {
			java.util.Date origin = new java.util.Date(0);
			
			double d = 0.0;
			
			try {
				int order = orderCombo.getSelectedIndex();
				boolean fromLast = false;
				if (order > 3) {
					fromLast = true;
					order = 8 - order - 1;
				}
			
				String prefix = "";
				if (prefixCheckBox.isSelected()) {
					prefix = prefixText.getText();
					d = guessDateFromPrefix(data.taxonList.getTaxonId(i), prefix, order, fromLast);
				} else {
					d = guessDateFromOrder(data.taxonList.getTaxonId(i), order, fromLast);
				}
				
			} catch (NumberFormatException nfe) {
			}	
			
			if (offsetCheck.isSelected()) {
				double offset = offsetText.getValue().doubleValue();
				if (unlessCheck.isSelected()) {
					double unless = unlessText.getValue().doubleValue();
					double offset2 = offset2Text.getValue().doubleValue();
					if (d < unless) {
						d += offset2;
					} else {
						d += offset;
					}
				} else {
					d += offset;
				}
			}	
				
			Date date = Date.createTimeSinceOrigin(d, Units.YEARS, origin);
			data.taxonList.getTaxon(i).setAttribute("date", date);
		}	
						
		// adjust the dates to the current timescale...	
		timeScaleChanged();

		dataTableModel.fireTableDataChanged();
		frame.fireDataChanged();
	}
	
	public double guessDateFromOrder(String label, int order, boolean fromLast) throws NumberFormatException {
		
		ArrayList fields = new ArrayList();
		
		int i = 0; 
		
		char c = label.charAt(i);
		
		do {
			// first find a part of a number
			while (!Character.isDigit(c) && c != '.') {
				i++;
				if (i == label.length()) break;
				c = label.charAt(i);
			}
			int j = i;
			
			if (i < label.length()) {
			
				// now find the end of the number
				while (Character.isDigit(c) || c == '.') {
					i++;
					if (i == label.length()) break;
					c = label.charAt(i);
				}
				
				fields.add(label.substring(j, i));
			}
									
		} while (i < label.length());
					
		int index = -1;
		
		if (fromLast) {
			index = fields.size() - order - 1;
		} else {
			index = order;
		}
		
		if (index < 0 || index >= fields.size()) {
			throw new NumberFormatException("Missing number field in taxon label");
		}
		
		return Double.parseDouble((String)fields.get(index));
	}
	
 	public double guessDateFromPrefix(String label, String prefix, int order, boolean fromLast) throws NumberFormatException {
		
		ArrayList fields = new ArrayList();
		
		String subLabel = label;
		
		int i = subLabel.indexOf(prefix);
		while (i != -1) {
			i += prefix.length();
			int j = i;
			
			if (i < label.length()) {
				// now find the beginning of the number
				char c = subLabel.charAt(i);
				while (i < subLabel.length() && (Character.isDigit(c) || c == '.')) {
					i++;
					if (i == subLabel.length()) break;
					c = subLabel.charAt(i);
				}

				fields.add(subLabel.substring(j, i));
			}
						
			subLabel = subLabel.substring(i);
			i = subLabel.indexOf(prefix);
		}
		
		int index = -1;
		
		if (fromLast) {
			index = fields.size() - order - 1;
		} else {
			index = order;
		}
		
		if (index < 0 || index >= fields.size()) {
			new NumberFormatException("Missing number field in taxon label");
		}
		
		return Double.parseDouble((String)fields.get(index));
	}
	
 	public class ClearDatesAction extends AbstractAction {
  		public ClearDatesAction() {
  			super("Clear Dates");
   			setToolTipText("Use this tool to add zero sampling dates to each taxon");
 		}
  		
  		public void actionPerformed(ActionEvent ae) { clearDates(); }
  	};

  	public class GuessDatesAction extends AbstractAction {
  		public GuessDatesAction() {
  			super("Guess Dates");
   			setToolTipText("Use this tool to guess the sampling dates from the taxon labels");
 		}
  		
  		public void actionPerformed(ActionEvent ae) { guessDates(); }
  	};

	class DataTableModel extends AbstractTableModel {

		String[] columnNames = { "Name", "Date" };

		public DataTableModel() {
		}

		public int getColumnCount() { 
			return columnNames.length; 
		}
		
		public int getRowCount() {
			return data.taxonList.getTaxonCount();
		}
		
		public Object getValueAt(int row, int col) {	
			switch (col) {
				case 0: return data.taxonList.getTaxonId(row); 
				case 1:
					Date date = data.taxonList.getTaxon(row).getDate();
					if (date != null) { 
						return new Double(date.getTimeValue());
					} else {
						return "-";
					} 
			}
			return null;
		}

		public void setValueAt(Object aValue, int row, int col) {	
			if (col == 0) {
				data.taxonList.getTaxon(row).setId(aValue.toString());
			} else if (col == 1) {
				Date date = data.taxonList.getTaxon(row).getDate();
				if (date != null) { 
					double d = ((Double)aValue).doubleValue();
					Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
					data.taxonList.getTaxon(row).setDate(newDate);
				} 
			}
			
			dataChanged();
		}

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return true;
            if (col == 1) {
				Date date = data.taxonList.getTaxon(row).getDate();
				return (date != null);
			}
			return false;
        }

		public String getColumnName(int column) {
			return columnNames[column]; 
		}
		
		public Class getColumnClass(int c) {return getValueAt(0, c).getClass();}
		
		public String toString() {
			StringBuffer buffer = new StringBuffer();

			buffer.append(getColumnName(0));
			for (int j = 1; j < getColumnCount(); j++) {
				buffer.append("\t");
				buffer.append(getColumnName(j));
			}
			buffer.append("\n");
			
			for (int i = 0; i < getRowCount(); i++) {
				buffer.append(getValueAt(i, 0));
				for (int j = 1; j < getColumnCount(); j++) {
					buffer.append("\t");
					buffer.append(getValueAt(i, j));
				}
				buffer.append("\n");
			}
			
			return buffer.toString();
		}
	};

}
