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

package dr.app.beauti;

import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.*;
import dr.evolution.util.*;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.*;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import java.awt.*;
import java.awt.event.*;
import java.util.Calendar;
import java.util.GregorianCalendar;

/**
 * @author			Andrew Rambaut
 * @author			Alexei Drummond
 * @version			$Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class DataPanel extends JPanel implements Exportable {

	JScrollPane scrollPane = new JScrollPane();
	JTable dataTable = null;
	DataTableModel dataTableModel = null;

	ClearDatesAction clearDatesAction = new ClearDatesAction();
	GuessDatesAction guessDatesAction = new GuessDatesAction();

	JComboBox unitsCombo = new JComboBox(new String[] {"Years", "Months", "Days"});
	JComboBox directionCombo = new JComboBox(new String[] {"Since some time in the past", "Before the present"});
	//RealNumberField originField = new RealNumberField(0.0, Double.POSITIVE_INFINITY);

	JComboBox translationCombo = new JComboBox();

	TableRenderer sequenceRenderer = null;

	BeautiFrame frame = null;

	BeautiOptions options = null;

	double[] heights = null;

	public DataPanel(BeautiFrame parent) {

		this.frame = parent;

		dataTableModel = new DataTableModel();
		TableSorter sorter = new TableSorter(dataTableModel);
		dataTable = new JTable(sorter);

		sorter.setTableHeader(dataTable.getTableHeader());

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

		dataTable.getColumnModel().getColumn(2).setCellRenderer(
				new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
		dataTable.getColumnModel().getColumn(2).setPreferredWidth(80);

		sequenceRenderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
		sequenceRenderer.setFont(new Font("Courier", Font.PLAIN, 12));
		dataTable.getColumnModel().getColumn(3).setCellRenderer(sequenceRenderer);
		dataTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);

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
		directionCombo.setOpaque(false);
		directionCombo.setEnabled(false);
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


		translationCombo.setOpaque(false);
		translationCombo.addItem("None");
		for (int i = 0; i < GeneticCode.GENETIC_CODE_DESCRIPTIONS.length; i++) {
			translationCombo.addItem(GeneticCode.GENETIC_CODE_DESCRIPTIONS[i]);
		}
		translationCombo.setEnabled(false);
		translationCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) { translationChanged(); }
		});

		JToolBar toolBar2 = new JToolBar();
		toolBar2.setOpaque(false);
		toolBar2.setFloatable(false);
		toolBar2.setLayout(new FlowLayout(java.awt.FlowLayout.LEFT, 0, 0));
		toolBar2.add(new JLabel("Translation:"));
		toolBar2.add(translationCombo);

		setOpaque(false);
		setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
		setLayout(new BorderLayout(0,0));
		add(toolBar1, "North");
		add(scrollPane, "Center");
		add(toolBar2, "South");

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
		calculateHeights();
		frame.dataChanged();
	}

	public final void timeScaleChanged() {
		int units = Units.YEARS;
		switch (unitsCombo.getSelectedIndex()) {
			case 0: units = Units.YEARS; break;
			case 1: units = Units.MONTHS; break;
			case 2: units = Units.DAYS; break;
		}

		boolean backwards = directionCombo.getSelectedIndex() == 1;

		//double origin = originField.getValue().doubleValue();

		for (int i = 0; i < options.alignment.getTaxonCount(); i++) {
			Date date = options.alignment.getTaxon(i).getDate();
			double d = date.getTimeValue();

			Date newDate = createDate(d, units, backwards, 0.0);

			options.alignment.getTaxon(i).setDate(newDate);
		}

		calculateHeights();

		dataTableModel.fireTableDataChanged();
		frame.dataChanged();
	}

	private Date createDate(double timeValue, int units, boolean backwards, double origin) {
		if (backwards) {
			return Date.createTimeAgoFromOrigin(timeValue, units, origin);
		} else {
			return Date.createTimeSinceOrigin(timeValue, units, origin);
		}
	}

	public final void translationChanged() {
		int index = translationCombo.getSelectedIndex() - 1;

		if (index < 0) {
			options.alignment = options.originalAlignment;
		} else {
			options.alignment = new ConvertAlignment(AminoAcids.INSTANCE, GeneticCode.GENETIC_CODES[index],
					options.originalAlignment);
		}

		dataTableModel.fireTableDataChanged();
		frame.dataChanged();
	}

	public void setOptions(BeautiOptions options) {

		this.options = options;

		if (options.originalAlignment != null) {
			clearDatesAction.setEnabled(true);
			guessDatesAction.setEnabled(true);
			unitsCombo.setEnabled(true);
			directionCombo.setEnabled(true);

			//originField.setEnabled(true);

			if (options.originalAlignment.getDataType() == Nucleotides.INSTANCE) {
				translationCombo.setEnabled(true);
				translationCombo.setSelectedIndex(options.translation);
			} else {
				translationCombo.setEnabled(false);
				translationCombo.setSelectedIndex(0);
			}

			sequenceRenderer.setText(options.alignment.getSequence(0).getSequenceString());
			int w = sequenceRenderer.getPreferredSize().width + 8;
			dataTable.getColumnModel().getColumn(3).setPreferredWidth(w);
		}

		unitsCombo.setSelectedIndex(options.datesUnits);
		directionCombo.setSelectedIndex(options.datesDirection);

		calculateHeights();

		dataTableModel.fireTableDataChanged();
	}

	public void getOptions(BeautiOptions options) {
		options.datesUnits = unitsCombo.getSelectedIndex();
		options.datesDirection = directionCombo.getSelectedIndex();
		options.translation = translationCombo.getSelectedIndex();
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
		int option = JOptionPane.showConfirmDialog(this, "Are you sure you wish to delete\n"+
				"the selected sequences?\n"+
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
				int index = options.originalAlignment.getTaxonIndex(names[i]);
				options.originalAlignment.removeSequence(index);
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

	}

	public void clearDates() {
		for (int i = 0; i < options.originalAlignment.getTaxonCount(); i++) {
			java.util.Date origin = new java.util.Date(0);

			double d = 0.0;

			Date date = Date.createTimeSinceOrigin(d, Units.YEARS, origin);
			options.originalAlignment.getTaxon(i).setAttribute("date", date);
		}

		// adjust the dates to the current timescale...
		timeScaleChanged();

		dataTableModel.fireTableDataChanged();
		frame.dataChanged();
	}

	public void guessDates() {

		OptionsPanel optionPanel = new OptionsPanel();

		optionPanel.addLabel("The date is given by a numerical field in the taxon label that is:");

		final JRadioButton orderRadio = new JRadioButton("Defined by its order", true);
		final JComboBox orderCombo = new JComboBox(new String[] {"first", "second", "third",
				"fourth", "fourth from last",
				"third from last", "second from last", "last"});

		optionPanel.addComponents(orderRadio, orderCombo);
		optionPanel.addSeparator();

		final JRadioButton prefixRadio = new JRadioButton("Defined by a prefix", false);
		final JTextField prefixText = new JTextField(16);
		prefixText.setEnabled(false);
		optionPanel.addComponents(prefixRadio, prefixText);
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
		Calendar calendar = GregorianCalendar.getInstance();

		int year = calendar.get(Calendar.YEAR) - 1999;
		unlessText.setValue(year);
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

		ButtonGroup group = new ButtonGroup();
		group.add(orderRadio);
		group.add(prefixRadio);
		ItemListener listener = new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				orderCombo.setEnabled(orderRadio.isSelected());
				prefixText.setEnabled(prefixRadio.isSelected());
			}
		};
		orderRadio.addItemListener(listener);
		prefixRadio.addItemListener(listener);

		JOptionPane optionPane = new JOptionPane(optionPanel,
				JOptionPane.QUESTION_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				null,
				null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		JDialog dialog = optionPane.createDialog(frame, "Guess Dates");
//		dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		dialog.setVisible(true);

		if (optionPane.getValue() == null) {
			return;
		}

		int value = ((Integer)optionPane.getValue()).intValue();
		if (value == -1 || value == JOptionPane.CANCEL_OPTION) {
			return;
		}

		for (int i = 0; i < options.originalAlignment.getTaxonCount(); i++) {
			java.util.Date origin = new java.util.Date(0);

			double d = 0.0;

			try {
				if (orderRadio.isSelected()) {
					int order = orderCombo.getSelectedIndex();
					boolean fromLast = false;
					if (order > 3) {
						fromLast = true;
						order = 8 - order - 1;
					}

					d = guessDateFromOrder(options.originalAlignment.getTaxonId(i), order, fromLast);
				} else {
					String prefix = prefixText.getText();
					d = guessDateFromPrefix(options.originalAlignment.getTaxonId(i), prefix);
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
			options.originalAlignment.getTaxon(i).setAttribute("date", date);
		}

		// adjust the dates to the current timescale...
		timeScaleChanged();

		dataTableModel.fireTableDataChanged();
		frame.dataChanged();
	}

	public double guessDateFromOrder(String label, int order, boolean fromLast) throws NumberFormatException {

		String field = null;

		if (fromLast) {
			int count = 0;
			int i = label.length() - 1;

			char c = label.charAt(i);

			do {
				// first find a part of a number
				while (!Character.isDigit(c) && c != '.') {
					i--;
					if (i < 0) break;
					c = label.charAt(i);
				}

				if (i < 0) new NumberFormatException("Missing number field in taxon label");

				int j = i + 1;

				// now find the beginning of the number
				while (Character.isDigit(c) || c == '.') {
					i--;
					if (i < 0) break;
					c = label.charAt(i);
				}

				field = label.substring(i+1, j);

				count++;

			} while (count <= order);

		} else {
			int count = 0;
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

				if (i == label.length()) new NumberFormatException("Missing number field in taxon label");

				// now find the beginning of the number
				while (Character.isDigit(c) || c == '.') {
					i++;
					if (i == label.length()) break;
					c = label.charAt(i);
				}

				field = label.substring(j, i);

				count++;

			} while (count <= order);
		}

		return Double.parseDouble(field);
	}

	public double guessDateFromPrefix(String label, String prefix) throws NumberFormatException {

		int i = label.indexOf(prefix);

		if (i == -1) new NumberFormatException("Missing prefix in taxon label");

		i += prefix.length();
		int j = i;

		// now find the beginning of the number
		char c = label.charAt(i);
		while (i < label.length() && (Character.isDigit(c) || c == '.')) {
			i++;
			c = label.charAt(i);
		}

		if (i == j) new NumberFormatException("Missing field after prefix in taxon label");

		String field = label.substring(j, i);

		return Double.parseDouble(field);
	}

	public class ClearDatesAction extends AbstractAction {
		public ClearDatesAction() {
			super("Clear Dates");
			setToolTipText("Use this tool to remove sampling dates from each taxon");
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

	private void calculateHeights() {

		options.maximumTipHeight = 0.0;
		if (options.alignment == null) return;

		heights = null;

		dr.evolution.util.Date mostRecent = null;
		for (int i = 0; i < options.alignment.getSequenceCount(); i++) {
			Date date = options.alignment.getTaxon(i).getDate();
			if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
				mostRecent = date;
			}
		}

		if (mostRecent != null) {
			heights = new double[options.alignment.getSequenceCount()];

			TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
			double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

			for (int i = 0; i < options.alignment.getSequenceCount(); i++) {
				Date date = options.alignment.getTaxon(i).getDate();
				if (date != null) {
					heights[i] = timeScale.convertTime(date.getTimeValue(), date) - time0;
					if (heights[i] > options.maximumTipHeight) options.maximumTipHeight = heights[i];
				}
			}
		}
	}

	class DataTableModel extends AbstractTableModel {

		String[] columnNames = { "Name", "Date", "Height", "Sequence" };

		public DataTableModel() {
		}

		public int getColumnCount() {
			return columnNames.length;
		}

		public int getRowCount() {
			if (options == null) return 0;
			if (options.alignment == null) return 0;

			return options.alignment.getTaxonCount();
		}

		public Object getValueAt(int row, int col) {
			switch (col) {
				case 0: return options.alignment.getTaxonId(row);
				case 1:
					Date date = options.alignment.getTaxon(row).getDate();
					if (date != null) {
						return new Double(date.getTimeValue());
					} else {
						return "-";
					}
				case 2:
					if (heights != null) {
						return new Double(heights[row]);
					} else {
						return "0.0";
					}
				case 3: return options.alignment.getAlignedSequenceString(row);
			}
			return null;
		}

		public void setValueAt(Object aValue, int row, int col) {
			if (col == 0) {
				options.alignment.getTaxon(row).setId(aValue.toString());
			} else if (col == 1) {
				Date date = options.alignment.getTaxon(row).getDate();
				if (date != null) {
					double d = ((Double)aValue).doubleValue();
					Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
					options.alignment.getTaxon(row).setDate(newDate);
				}
			}

			dataChanged();
		}

		public boolean isCellEditable(int row, int col) {
			if (col == 0) return true;
			if (col == 1) {
				Date date = options.alignment.getTaxon(row).getDate();
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
