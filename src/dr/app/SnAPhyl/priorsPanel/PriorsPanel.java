/*
 * PriorsPanel.java
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

package dr.app.SnAPhyl.priorsPanel;

import dr.app.SnAPhyl.BeautiFrame;

import dr.app.beauti.BeautiPanel;
import dr.app.beauti.enumTypes.TreePriorType;
import dr.app.beauti.options.*;
import dr.app.beauti.priorsPanel.DiscretePriorDialog;
import dr.app.beauti.priorsPanel.PriorDialog;

import dr.util.NumberFormatter;
import org.virion.jam.framework.Exportable;
import org.virion.jam.panels.OptionsPanel;
import org.virion.jam.table.HeaderRenderer;
import org.virion.jam.table.TableEditorStopper;
import org.virion.jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.EnumSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorsPanel extends BeautiPanel implements Exportable {
    
    private static final long serialVersionUID = -2936049032365493416L;
    JScrollPane scrollPane = new JScrollPane();
    JTable priorTable = null;
    PriorTableModel priorTableModel = null;

    public ArrayList<Parameter> parameters = new ArrayList<Parameter>();

    BeautiFrame frame = null;
    BeautiOptions options = null;

    private final boolean isDefaultOnly;

    private JComboBox treePriorCombo = new JComboBox(EnumSet.range(TreePriorType.SPECIES_YULE, TreePriorType.SPECIES_BIRTH_DEATH).toArray());

    public PriorsPanel(BeautiFrame parent, boolean isDefaultOnly) {
        this.frame = parent;
        this.isDefaultOnly = isDefaultOnly;

        priorTableModel = new PriorTableModel(this);
        priorTable = new JTable(priorTableModel);

        priorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        priorTable.getTableHeader().setReorderingAllowed(false);
        priorTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        priorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_LAST_COLUMN);
        priorTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(0).setPreferredWidth(160);

        priorTable.getColumnModel().getColumn(1).setCellRenderer(
                new ButtonRenderer(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
        priorTable.getColumnModel().getColumn(1).setCellEditor(
                new ButtonEditor(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
        priorTable.getColumnModel().getColumn(1).setPreferredWidth(260);

        priorTable.getColumnModel().getColumn(2).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        if (isDefaultOnly) {
            priorTable.getColumnModel().getColumn(2).setPreferredWidth(60);
        } else {
            priorTable.getColumnModel().getColumn(2).setPreferredWidth(30);
        }

        priorTable.getColumnModel().getColumn(3).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(3).setPreferredWidth(400);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(priorTable);

        scrollPane = new JScrollPane(priorTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

        scrollPane.setOpaque(false);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.setOpaque(false);

        if (isDefaultOnly) {
            scrollPane.setPreferredSize(new java.awt.Dimension(800, 600));
            panel.add(scrollPane, BorderLayout.CENTER);
            panel.add(new JLabel("These priors listed above are still set to the default values " +
                    "and need to be reviewed, especially their upper and lower limits."), BorderLayout.SOUTH);
        } else {
            panel.add(new JLabel("Priors for model parameters and statistics:"), BorderLayout.NORTH);
            panel.add(scrollPane, BorderLayout.CENTER);
            panel.add(new JLabel("* Marked parameters currently have a default prior distribution. " +
                    "You should check that these are appropriate."), BorderLayout.SOUTH);

            OptionsPanel treePriorPanel = new OptionsPanel(12, 18);
            treePriorPanel.addComponentWithLabel("Species Tree Prior:", treePriorCombo);
            treePriorCombo.addItemListener(
                new ItemListener() {
                    public void itemStateChanged(ItemEvent ev) {
//                    	partitionTreePrior.setNodeHeightPrior( (TreePriorType) treePriorCombo.getSelectedItem());
                    }
                }
            );

            add(treePriorPanel, BorderLayout.NORTH);
                        
        }

        add(panel, BorderLayout.CENTER);
    }

    public void setOptions(BeautiOptions options) {
    	this.options = options;
    	
        parameters = options.selectParameters();
        priorTableModel.fireTableDataChanged();

        validate();
        repaint();
    }

    public void setParametersList(BeautiOptions options) {
        this.options = options;

        parameters.clear();
        for (Parameter param : options.selectParameters()) {
            if (!param.isPriorEdited()) {
                parameters.add(param);
            }
        }
    }

    private PriorDialog priorDialog = null;
    private DiscretePriorDialog discretePriorDialog = null;

    private void priorButtonPressed(int row) {
        Parameter param = parameters.get(row);

        if (param.isDiscrete) {
            if (discretePriorDialog == null) {
                discretePriorDialog = new DiscretePriorDialog(frame);
            }

            if (discretePriorDialog.showDialog(param) == JOptionPane.CANCEL_OPTION) {
                return;
            }
        } else {
            if (priorDialog == null) {
                priorDialog = new PriorDialog(frame);
            }

            if (priorDialog.showDialog(param) == JOptionPane.CANCEL_OPTION) {
                return;
            }
        }
        param.setPriorEdited(true);

        if (isDefaultOnly) {             
            setParametersList(options);
        }
        
        if (param.getBaseName().endsWith("treeModel.rootHeight") || param.taxa != null) { // param.taxa != null is TMRCA 
        	if (options.clockModelOptions.isNodeCalibrated(param)) {
        		options.clockModelOptions.nodeCalibration();
        		frame.setAllOptions();
        		frame.setStatusMessage();
//        	} else {
//        		options.clockModelOptions.fixRateOfFirstClockPartition();
        	}
        }

        priorTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
    }

    public JComponent getExportableComponent() {
        return priorTable;
    }

    NumberFormatter formatter = new NumberFormatter(4);

    class DoubleRenderer extends TableRenderer {

        private static final long serialVersionUID = -2614341608257369805L;

        public DoubleRenderer(int alignment, Insets insets) {

            super(true, alignment, insets);
        }

        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected,
                                                       boolean hasFocus, int row, int column) {

            String s;
            if (((Double) value).isNaN()) {
                s = "random";
            } else {
                s = formatter.format((Double) value);
            }
            return super.getTableCellRendererComponent(table, s, isSelected, hasFocus, row, column);

        }
    }

    public class ButtonRenderer extends JButton implements TableCellRenderer {

        private static final long serialVersionUID = -2416184092883649169L;

        public ButtonRenderer(int alignment, Insets insets) {
            setOpaque(true);
            setHorizontalAlignment(alignment);
            setMargin(insets);
        }

        public Component getTableCellRendererComponent(JTable table, Object value,
                                                       boolean isSelected, boolean hasFocus, int row, int column) {
            setEnabled(table.isEnabled());
            setFont(table.getFont());
            if (isSelected) {
                setForeground(table.getSelectionForeground());
                setBackground(table.getSelectionBackground());
            } else {
                setForeground(table.getForeground());
                setBackground(UIManager.getColor("Button.background"));
            }
            setText((value == null) ? "" : value.toString());
            return this;
        }
    }

    public class ButtonEditor extends DefaultCellEditor {
        
        private static final long serialVersionUID = 6372738480075411674L;
        protected JButton button;
        private String label;
        private boolean isPushed;
        private int row;

        public ButtonEditor(int alignment, Insets insets) {
            super(new JCheckBox());
            button = new JButton();
            button.setOpaque(true);
            button.setHorizontalAlignment(alignment);
            button.setMargin(insets);
            button.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    fireEditingStopped();
                }
            });
        }

        public Component getTableCellEditorComponent(JTable table, Object value,
                                                     boolean isSelected, int row, int column) {
            button.setEnabled(table.isEnabled());
            button.setFont(table.getFont());
            if (isSelected) {
                button.setForeground(table.getSelectionForeground());
                button.setBackground(table.getSelectionBackground());
            } else {
                button.setForeground(table.getForeground());
                button.setBackground(table.getBackground());
            }
            label = (value == null) ? "" : value.toString();
            button.setText(label);
            isPushed = true;
            this.row = row;
            return button;
        }

        public Object getCellEditorValue() {
            if (isPushed) {
                priorButtonPressed(row);
            }
            isPushed = false;
            return label;
        }

        public boolean stopCellEditing() {
            isPushed = false;
            return super.stopCellEditing();
        }

        protected void fireEditingStopped() {
            super.fireEditingStopped();
        }
    }
}
