/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.SnAPhyl.datapanel;

import dr.app.SnAPhyl.SnAPhylFrame;

import dr.app.beauti.BeautiPanel;
import dr.app.beauti.alignmentviewer.*;
import dr.app.beauti.datapanel.BeautiAlignmentBuffer;
import dr.app.beauti.options.*;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.util.PanelUtils;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.DataType;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.*;
import dr.gui.table.DateCellEditor;
import dr.gui.table.TableSorter;
import org.virion.jam.components.RealNumberField;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.ActionPanel;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.Date;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DataPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class DataPanel extends BeautiPanel implements Exportable {

        JScrollPane scrollPane = new JScrollPane();
    JTable dataTable = null;
    DataTableModel dataTableModel = null;

    JComboBox unitsCombo = new JComboBox(new String[] {"Years", "Months", "Days"});
    JComboBox directionCombo = new JComboBox(new String[] {"Since some time in the past", "Before the present"});
    //RealNumberField originField = new RealNumberField(0.0, Double.POSITIVE_INFINITY);


    TableRenderer sequenceRenderer = null;

    SnAPhylFrame frame = null;

    dr.app.oldbeauti.BeautiOptions options = null;

    double[] heights = null;

    public DataPanel(SnAPhylFrame parent) {

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

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(dataTable);

        sequenceRenderer = new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4));
        sequenceRenderer.setFont(new Font("Courier", Font.PLAIN, 12));

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) { selectionChanged(); }
        });

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);


	    setupComponent(unitsCombo);
        unitsCombo.setEnabled(false);
	    setupComponent(directionCombo);
        directionCombo.setEnabled(false);
        //originField.setEnabled(false);
        //originField.setValue(0.0);
        //originField.setColumns(12);





        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0,0));

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

	private void setupComponent(JComponent comp) {
		comp.setOpaque(false);

		//comp.setFont(UIManager.getFont("SmallSystemFont"));
		//comp.putClientProperty("JComponent.sizeVariant", "small");
		if (comp instanceof JButton) {
			comp.putClientProperty("JButton.buttonType", "roundRect");
		}
		if (comp instanceof JComboBox) {
			comp.putClientProperty("JComboBox.isSquare", Boolean.TRUE);
		}

	}



    public final void timeScaleChanged() {
        Units.Type units = Units.Type.YEARS;
        switch (unitsCombo.getSelectedIndex()) {
            case 0: units = Units.Type.YEARS; break;
            case 1: units = Units.Type.MONTHS; break;
            case 2: units = Units.Type.DAYS; break;
        }

        boolean backwards = directionCombo.getSelectedIndex() == 1;

        //double origin = originField.getValue().doubleValue();

        for (int i = 0; i < options.taxonList.getTaxonCount(); i++) {
            dr.evolution.util.Date date = options.taxonList.getTaxon(i).getDate();
            double d = date.getTimeValue();

            dr.evolution.util.Date newDate = createDate(d, units, backwards, 0.0);

            options.taxonList.getTaxon(i).setDate(newDate);
        }



        dataTableModel.fireTableDataChanged();

    }

    private dr.evolution.util.Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return dr.evolution.util.Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return dr.evolution.util.Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }



    public void setOptions(dr.app.oldbeauti.BeautiOptions options) {

        this.options = options;

        if (options.taxonList != null) {

            unitsCombo.setEnabled(true);
            directionCombo.setEnabled(true);

            //originField.setEnabled(true);


        }

        setupTable();

        unitsCombo.setSelectedIndex(options.datesUnits);
        directionCombo.setSelectedIndex(options.datesDirection);



        dataTableModel.fireTableDataChanged();
    }

    private void setupTable() {

        dataTableModel.fireTableStructureChanged();
        if (options.alignment != null) {

            dataTable.getColumnModel().getColumn(3).setCellRenderer(sequenceRenderer);
            dataTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);

            sequenceRenderer.setText(options.alignment.getSequence(0).getSequenceString());
            int w = sequenceRenderer.getPreferredSize().width + 8;
            dataTable.getColumnModel().getColumn(3).setPreferredWidth(w);
        }
    }

    public void getOptions(dr.app.oldbeauti.BeautiOptions options) {
        options.datesUnits = unitsCombo.getSelectedIndex();
        options.datesDirection = directionCombo.getSelectedIndex();

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
                if (options.originalAlignment != null) {
                    int index = options.originalAlignment.getTaxonIndex(names[i]);
                    options.originalAlignment.removeSequence(index);
                } else {
                    // there is no alignment so options.taxonList must be a Taxa object:
                    int index = options.taxonList.getTaxonIndex(names[i]);
                    ((Taxa)options.taxonList).removeTaxon(options.taxonList.getTaxon(index));
                }
            }

            if (options.taxonList.getTaxonCount() == 0) {
                // if all the sequences are deleted we may as well throw
                // away the alignment...

                options.originalAlignment = null;
                options.alignment = null;
                options.taxonList = null;
            }

            dataTableModel.fireTableDataChanged();

        }

    }

   


    @Override
    public void setOptions(BeautiOptions options) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void getOptions(BeautiOptions options) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


    class DataTableModel extends AbstractTableModel {

        /**
         *
         */

        String[] columnNames1 = { "Name", "Date", "Height", "Sequence" };
        String[] columnNames2 = { "Name", "Date", "Height" };

        public DataTableModel() {
        }

        public int getColumnCount() {
            if (options != null && options.alignment != null) {
                return columnNames1.length;
            } else {
                return columnNames2.length;
            }
        }

        public int getRowCount() {
            if (options == null) return 0;
            if (options.taxonList == null) return 0;

            return options.taxonList.getTaxonCount();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0: return options.taxonList.getTaxonId(row);
                case 1:
                    dr.evolution.util.Date date = options.taxonList.getTaxon(row).getDate();
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
                options.taxonList.getTaxon(row).setId(aValue.toString());
            } else if (col == 1) {
                dr.evolution.util.Date date = options.taxonList.getTaxon(row).getDate();
                if (date != null) {
                    double d = ((Double)aValue).doubleValue();
                    dr.evolution.util.Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
                    options.taxonList.getTaxon(row).setDate(newDate);
                }
            }


        }

        public boolean isCellEditable(int row, int col) {
            if (col == 0) return true;
            if (col == 1) {
                dr.evolution.util.Date date = options.taxonList.getTaxon(row).getDate();
                return (date != null);
            }
            return false;
        }

        public String getColumnName(int column) {
            return columnNames1[column];
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