/*
 * MeasurementsPanel.java
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

package dr.app.mapper.application;

import dr.app.gui.FileDrop;
import dr.app.gui.table.DateCellEditor;
import dr.app.gui.table.TableEditorStopper;
import dr.app.gui.table.TableSorter;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.util.DataTable;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id: StrainsPanel.java,v 1.17 2006/09/05 13:29:34 rambaut Exp $
 */
public class MeasurementsPanel extends JPanel implements Exportable, MapperDocument.Listener {

    private JScrollPane scrollPane = new JScrollPane();
    private JTable dataTable = null;
    private DataTableModel dataTableModel = null;

    private final MapperFrame frame;
    private final MapperDocument document;


    public MeasurementsPanel(final MapperFrame parent, final MapperDocument document) {

        this.frame = parent;
        this.document = document;

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

        dataTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                selectionChanged();
            }
        });

        scrollPane = new JScrollPane(dataTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setOpaque(false);

        JToolBar toolBar1 = new JToolBar();
        toolBar1.setFloatable(false);
        toolBar1.setOpaque(false);

//        toolBar1.setLayout(new FlowLayout(FlowLayout.LEFT, 0, 0));
//        JButton button = new JButton(clearDatesAction);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);
//        button = new JButton(guessDatesAction);
//        PanelUtils.setupComponent(button);
//        toolBar1.add(button);
//        toolBar1.add(new JToolBar.Separator(new Dimension(12, 12)));
//        final JLabel unitsLabel = new JLabel("Dates specified as ");
//        toolBar1.add(unitsLabel);
//        toolBar1.add(unitsCombo);
//        toolBar1.add(directionCombo);

        setOpaque(false);
        setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        setLayout(new BorderLayout(0, 0));

        add(toolBar1, "North");
        add(scrollPane, "Center");

        Color focusColor = UIManager.getColor("Focus.color");
        Border focusBorder = BorderFactory.createMatteBorder(2, 2, 2, 2, focusColor);
        scrollPane.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
        new FileDrop(null, scrollPane, focusBorder, new FileDrop.Listener() {
            public void filesDropped(java.io.File[] files) {
                frame.importMeasurementFiles(files);
            }   // end filesDropped
        }); // end FileDrop.Listener

    }

    @Override
    public void taxaChanged() {
        dataTableModel.fireTableDataChanged();
    }

    public JComponent getExportableComponent() {
        return dataTable;
    }

    public void selectionChanged() {
        // nothing to do
    }

    class DataTableModel extends AbstractTableModel {

        String[] columnNames = {"Serum", "Virus", "Titre", "Table"};

        public DataTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            java.util.List<MapperDocument.Measurement> measurementList = document.getMeasurements();
            return measurementList.size();
        }

        public Object getValueAt(int row, int col) {
            MapperDocument.Measurement measurement = document.getMeasurements().get(row);

            switch (col) {
                case 0:
                    return measurement.columnStrain.getId();
                case 1:
                    return measurement.rowStrain.getId();
                case 2:
//                    if (measurement.type == MapperDocument.MeasurementType.THRESHOLD) {
//                        return "<" + measurement.titre;
//                    }
                    return measurement.titre;
                case 3:
                    return measurement.column;
            }
            return null;
        }

        public void setValueAt(Object aValue, int row, int col) {
//            java.util.List<Taxon> taxonList = document.getTaxa();
//
//            if (col == 0) {
//                taxonList.get(row).setId(aValue.toString());
//            } else if (col == 1) {
//                Date date = taxonList.get(row).getDate();
//                if (date != null) {
//                    double d = (Double) aValue;
//                    Date newDate = createDate(d, date.getUnits(), date.isBackwards(), date.getOrigin());
//                    taxonList.get(row).setDate(newDate);
//                }
//            }
//
//            timeScaleChanged();
        }

        public boolean isCellEditable(int row, int col) {
//            if (col == 0) return true;
//            if (col == 1) {
//                Date date = document.getTaxa().get(row).getDate();
//                return (date != null);
//            }
            return false;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

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
    }
}
