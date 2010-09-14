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

package dr.app.beauti.priorsPanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.BeautiPanel;
import dr.app.beauti.enumTypes.ClockType;
import dr.app.beauti.enumTypes.PriorType;
import dr.app.beauti.options.*;
import dr.app.gui.table.TableEditorStopper;
import dr.util.NumberFormatter;
import jam.framework.Exportable;
import jam.table.HeaderRenderer;
import jam.table.TableRenderer;

import javax.swing.*;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;

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
    JLabel messageLabel = new JLabel();

    public ArrayList<Parameter> parameters = new ArrayList<Parameter>();

    BeautiFrame frame = null;
    BeautiOptions options = null;

    boolean hasUndefinedPrior = false;
    boolean hasImproperPrior = false;
    boolean hasRate = false;

    private final boolean isDefaultOnly;

    public PriorsPanel(BeautiFrame parent, boolean isDefaultOnly) {
        this.frame = parent;
        this.isDefaultOnly = isDefaultOnly;

        priorTableModel = new PriorTableModel(this);
        priorTable = new JTable(priorTableModel);

        priorTable.setAutoResizeMode(javax.swing.JTable.AUTO_RESIZE_OFF);
        priorTable.getTableHeader().setReorderingAllowed(false);
        priorTable.getTableHeader().setDefaultRenderer(
                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        priorTable.getColumnModel().getColumn(0).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(0).setMinWidth(200);

        priorTable.getColumnModel().getColumn(1).setCellRenderer(
                new ButtonRenderer(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
        priorTable.getColumnModel().getColumn(1).setCellEditor(
                new ButtonEditor(SwingConstants.LEFT, new Insets(0, 8, 0, 8)));
        priorTable.getColumnModel().getColumn(1).setMinWidth(260);

        priorTable.getColumnModel().getColumn(2).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(2).setMinWidth(30);


        priorTable.getColumnModel().getColumn(3).setCellRenderer(
                new TableRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));
        priorTable.getColumnModel().getColumn(3).setMinWidth(410);

        TableEditorStopper.ensureEditingStopWhenTableLosesFocus(priorTable);

        scrollPane = new JScrollPane(priorTable,
                JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

        scrollPane.setOpaque(false);

        messageLabel.setText(getMessage());

        setOpaque(false);
        setLayout(new BorderLayout(0, 0));
        setBorder(new BorderUIResource.EmptyBorderUIResource(new java.awt.Insets(12, 12, 12, 12)));

        if (!isDefaultOnly) {
            add(new JLabel("Priors for model parameters and statistics:"), BorderLayout.NORTH);
        }

        add(scrollPane, BorderLayout.CENTER);
        add(messageLabel, BorderLayout.SOUTH);
    }

    public void setOptions(BeautiOptions options) {
        this.options = options;

        parameters = options.selectParameters();

        messageLabel.setText(getMessage());

        priorTableModel.fireTableDataChanged();

        validate();
        repaint();
    }

    private String getMessage() {
        String message = "<html>";
        if (isDefaultOnly) {
            message += "<ul><li>These priors listed above are still set to the default values " +
                    "and need to be reviewed.</li>";

            hasUndefinedPrior = false;
            hasImproperPrior = false;
            hasRate = false;

            for (Parameter param : parameters) {
                if (param.priorType == PriorType.UNDEFINED) {
                    hasUndefinedPrior = true;
                }
                if (param.isPriorImproper()) {
                    hasImproperPrior = true;
                }
                if (param.getBaseName().endsWith("clock.rate") || param.getBaseName().endsWith(ClockType.UCED_MEAN)
                        || param.getBaseName().endsWith(ClockType.UCLD_MEAN)) {
                    hasRate = true;
                }
            }

            if (hasUndefinedPrior) {
                message += "<li><b><font color=\"#E42217\">These priors need to be defined by the user.</font></b></li>";
            }
            if (hasImproperPrior) {
                message += "<li><b><font color=\"#B4B417\">Warning: one or more parameters have improper priors.</font></b></li>";
            }
            if (hasRate) {
                message += "<li>" +
                        "Priors on clock rates in particular should be proper such as a high variance lognormal or gamma distribution " +
                        "with a mean appropriate for the organism and units of time employed.</li>";
            }
            message += "</ul>";

        } else {
            message += "* Marked parameters currently have a default prior distribution. " +
                    "You should check that these are appropriate.";
        }


        return message + "</html>";
    }

    public void setParametersList(BeautiOptions options) {
        this.options = options;

        parameters.clear();
        for (Parameter param : options.selectParameters()) {
            if (!param.isPriorEdited() || param.isPriorImproper()) {
                parameters.add(param);
            }
        }
        messageLabel.setText(getMessage());
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
        protected Color undefinedColour = new Color(0xE4, 0x22, 0x17);
        protected Color improperColour = new Color(0xB4, 0xB4, 0x17);


        private static final long serialVersionUID = -2416184092883649169L;

        public ButtonRenderer(int alignment, Insets insets) {
            setOpaque(true);
            setHorizontalAlignment(alignment);
            setMargin(insets);

//            setFont(UIManager.getFont("SmallSystemFont"));
//            putClientProperty("JComponent.sizeVariant", "small");
//            putClientProperty("JButton.buttonType", "square");
        }

        public Component getTableCellRendererComponent(JTable table, Object value,
                                                       boolean isSelected, boolean hasFocus, int row, int column) {
            setEnabled(table.isEnabled());
            setFont(table.getFont());
            if (isSelected) {
                //setForeground(table.getSelectionForeground());
                setBackground(table.getSelectionBackground());
            } else {
                //setForeground(table.getForeground());
                setBackground(UIManager.getColor("Button.background"));
            }

            if (value.toString().startsWith("?")) {
                setForeground(undefinedColour);
            } else if (value.toString().startsWith("!")) {
                setForeground(improperColour);
            } else {
                setForeground(UIManager.getColor("Button.foreground"));
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
//            button.setFont(UIManager.getFont("SmallSystemFont"));
//            button.putClientProperty("JComponent.sizeVariant", "small");
//            button.putClientProperty("JButton.buttonType", "square");
        }

        public Component getTableCellEditorComponent(JTable table, Object value,
                                                     boolean isSelected, int row, int column) {
            button.setEnabled(table.isEnabled());
            button.setFont(table.getFont());
            if (isSelected) {
//                button.setForeground(table.getSelectionForeground());
                button.setBackground(table.getSelectionBackground());
            } else {
//                button.setForeground(table.getForeground());
                button.setBackground(UIManager.getColor("Button.background"));
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
