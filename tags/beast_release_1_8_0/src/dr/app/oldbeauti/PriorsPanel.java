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

package dr.app.oldbeauti;

import dr.app.gui.components.RealNumberField;
import dr.app.gui.table.TableEditorStopper;
import dr.util.NumberFormatter;
import dr.app.gui.components.WholeNumberField;
import jam.framework.Exportable;
import jam.panels.OptionsPanel;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: PriorsPanel.java,v 1.9 2006/09/05 13:29:34 rambaut Exp $
 */
public class PriorsPanel extends JPanel implements Exportable {

    /**
     *
     */
    private static final long serialVersionUID = -2936049032365493416L;
    JScrollPane scrollPane = new JScrollPane();
    JTable priorTable = null;
    PriorTableModel priorTableModel = null;

    OptionsPanel treePriorPanel = new OptionsPanel();
    JComboBox treePriorCombo;
    JComboBox parameterizationCombo = new JComboBox(new String[]{
            "Growth Rate", "Doubling Time"});
    JComboBox bayesianSkylineCombo = new JComboBox(new String[]{
            "Piecewise-constant", "Piecewise-linear"});
    WholeNumberField groupCountField = new WholeNumberField(2, Integer.MAX_VALUE);

    RealNumberField samplingProportionField = new RealNumberField(Double.MIN_VALUE, 1.0);

    JCheckBox upgmaStartingTreeCheck = new JCheckBox("Use UPGMA to construct a starting tree");

    public ArrayList parameters = new ArrayList();

    BeautiFrame frame = null;

    public PriorsPanel(BeautiFrame parent) {

        this.frame = parent;

        priorTableModel = new PriorTableModel();
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
        priorTable.getColumnModel().getColumn(2).setPreferredWidth(400);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(priorTable);

        scrollPane = new JScrollPane(priorTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

        scrollPane.setOpaque(false);

        java.awt.event.ItemListener listener = new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent ev) {
                if (!settingOptions) frame.priorsChanged();
            }
        };

        // order here must match corrosponding BeautiOptions constant, i.e. BeautiOptions.CONSTANT == 0 etc
        if (BeautiApp.developer) {
            treePriorCombo = new JComboBox(new String[]{
                    "Coalescent: Constant Size",
                    "Coalescent: Exponential Growth",
                    "Coalescent: Logistic Growth",
                    "Coalescent: Expansion Growth",
                    "Coalescent: Bayesian Skyline",
                    "Coalescent: Extended Bayesian Skyline",
                    "Speciation: Yule Process",
                    "Speciation: Birth-Death Process"
            });
        } else {
            treePriorCombo = new JComboBox(new String[]{
                    "Coalescent: Constant Size",
                    "Coalescent: Exponential Growth",
                    "Coalescent: Logistic Growth",
                    "Coalescent: Expansion Growth",
                    "Coalescent: Bayesian Skyline",
                    "Coalescent: Extended Bayesian Skyline",
                    "Speciation: Yule Process",
                    "Speciation: Birth-Death Process"
            });
        }
        setupComponent(treePriorCombo);
        treePriorCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        if (!settingOptions) frame.priorsChanged();
                        setupPanel();
                    }
                }
        );

        KeyListener keyListener = new KeyAdapter() {
            public void keyTyped(java.awt.event.KeyEvent ev) {
                if (!settingOptions && ev.getKeyCode() == KeyEvent.VK_ENTER) {
                    frame.priorsChanged();
                }
            }
        };

        groupCountField.addKeyListener(keyListener);
        samplingProportionField.addKeyListener(keyListener);

        FocusListener focusListener = new FocusAdapter() {
            public void focusLost(FocusEvent focusEvent) {
                frame.priorsChanged();
            }
        };
        samplingProportionField.addFocusListener(focusListener);
        groupCountField.addFocusListener(focusListener);

        setupComponent(parameterizationCombo);
        parameterizationCombo.addItemListener(listener);

        setupComponent(bayesianSkylineCombo);
        bayesianSkylineCombo.addItemListener(listener);

        setupComponent(upgmaStartingTreeCheck);

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        JPanel panel = new JPanel(new BorderLayout(0, 0));
        panel.setOpaque(false);
        panel.add(new JLabel("Priors for model parameters and statistics:"), BorderLayout.NORTH);
        panel.add(scrollPane, BorderLayout.CENTER);
        panel.add(new JLabel("* Marked parameters currently have a default prior distribution. " +
                "You should check that these are appropriate."), BorderLayout.SOUTH);

        treePriorPanel.setBorder(null);
        add(treePriorPanel, BorderLayout.NORTH);
        add(panel, BorderLayout.CENTER);
    }

    private void setupComponent(JComponent comp) {
        comp.setOpaque(false);

        if (comp instanceof JButton) {
            comp.putClientProperty("JButton.buttonType", "roundRect");
        }
        if (comp instanceof JComboBox) {
            comp.putClientProperty("JComboBox.isSquare", Boolean.TRUE);
        }
    }

    private void setupPanel() {

        treePriorPanel.removeAll();

        treePriorPanel.addComponentWithLabel("Tree Prior:", treePriorCombo);
        if (treePriorCombo.getSelectedIndex() == 1 || // exponential
                treePriorCombo.getSelectedIndex() == 2 || // logistic
                treePriorCombo.getSelectedIndex() == 3) { // expansion
            treePriorPanel.addComponentWithLabel("Parameterization for growth:", parameterizationCombo);
        } else if (treePriorCombo.getSelectedIndex() == 4) { // bayesian skyline
            groupCountField.setColumns(6);
            treePriorPanel.addComponentWithLabel("Number of groups:", groupCountField);
            treePriorPanel.addComponentWithLabel("Skyline Model:", bayesianSkylineCombo);
        } else if (treePriorCombo.getSelectedIndex() == 6) { // birth-death
            samplingProportionField.setColumns(8);
            treePriorPanel.addComponentWithLabel("Proportion of taxa sampled:", samplingProportionField);
        }

        treePriorPanel.addComponent(upgmaStartingTreeCheck);

        validate();
        repaint();
    }

    private boolean settingOptions = false;

    public void setOptions(BeautiOptions options) {
        settingOptions = true;
        parameters = options.selectParameters();
        priorTableModel.fireTableDataChanged();

        treePriorCombo.setSelectedIndex(options.nodeHeightPrior);

        groupCountField.setValue(options.skylineGroupCount);
        samplingProportionField.setValue(options.birthDeathSamplingProportion);

        parameterizationCombo.setSelectedIndex(options.parameterization);
        bayesianSkylineCombo.setSelectedIndex(options.skylineModel);

        upgmaStartingTreeCheck.setSelected(options.upgmaStartingTree);

        setupPanel();

        settingOptions = false;

        validate();
        repaint();
    }

    private PriorDialog priorDialog = null;
    private DiscretePriorDialog discretePriorDialog = null;

    private void priorButtonPressed(int row) {
        BeautiOptions.Parameter param = (BeautiOptions.Parameter) parameters.get(row);

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
        param.priorEdited = true;

        priorTableModel.fireTableDataChanged();
    }

    public void getOptions(BeautiOptions options) {
        if (settingOptions) return;

        if (treePriorCombo.getSelectedIndex() == BeautiOptions.CONSTANT) {
            options.nodeHeightPrior = BeautiOptions.CONSTANT;
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.EXPONENTIAL) {
            options.nodeHeightPrior = BeautiOptions.EXPONENTIAL;
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.LOGISTIC) {
            options.nodeHeightPrior = BeautiOptions.LOGISTIC;
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.EXPANSION) {
            options.nodeHeightPrior = BeautiOptions.EXPANSION;
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.SKYLINE) {
            options.nodeHeightPrior = BeautiOptions.SKYLINE;
            Integer groupCount = groupCountField.getValue();
            if (groupCount != null) {
                options.skylineGroupCount = groupCount;
            } else {
                options.skylineGroupCount = 5;
            }
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.EXTENDED_SKYLINE) {
            options.nodeHeightPrior = BeautiOptions.EXTENDED_SKYLINE;
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.YULE) {
            options.nodeHeightPrior = BeautiOptions.YULE;
        } else if (treePriorCombo.getSelectedIndex() == BeautiOptions.BIRTH_DEATH) {
            options.nodeHeightPrior = BeautiOptions.BIRTH_DEATH;
            Double samplingProportion = samplingProportionField.getValue();
            if (samplingProportion != null) {
                options.birthDeathSamplingProportion = samplingProportion;
            } else {
                options.birthDeathSamplingProportion = 1.0;
            }
        } else {
            throw new RuntimeException("Unexpected value from treePriorCombo");
        }

        options.parameterization = parameterizationCombo.getSelectedIndex();
        options.skylineModel = bayesianSkylineCombo.getSelectedIndex();

        options.upgmaStartingTree = upgmaStartingTreeCheck.isSelected();
    }

    public JComponent getExportableComponent() {
        return priorTable;
    }

    NumberFormatter formatter = new NumberFormatter(4);

    class PriorTableModel extends AbstractTableModel {

        /**
         *
         */
        private static final long serialVersionUID = -8864178122484971872L;
        String[] columnNames = {"Parameter", "Prior", "Description"};

        public PriorTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            return parameters.size();
        }

        public Object getValueAt(int row, int col) {
            BeastGenerator.Parameter param = (BeastGenerator.Parameter) parameters.get(row);
            switch (col) {
                case 0:
                    return param.getName();
                case 1:
                    return param.priorType.getPriorString(param);
                case 2:
                    return param.getDescription();
            }
            return null;
        }

        public String getColumnName(int column) {
            return columnNames[column];
        }

        public Class getColumnClass(int c) {
            return getValueAt(0, c).getClass();
        }

        public boolean isCellEditable(int row, int col) {
            return col == 1;
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

    class DoubleRenderer extends TableRenderer {

        /**
         *
         */
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

        /**
         *
         */
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
        /**
         *
         */
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
