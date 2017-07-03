/*
 * JointPriorDialog.java
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

/*
 * PriorDialog.java
 *
 * @author Marc A. Suchard
 */

package dr.app.beauti.priorspanel;

import dr.app.beauti.components.hpm.HierarchicalModelComponentOptions;
import dr.app.beauti.components.linkedparameters.LinkedParameter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.util.OSType;
import jam.panels.ActionPanel;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.util.List;

/**
 * A dialog which acts as a base for linking parameters together under joint priors such
 * as hierarchical models.
 *
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class JointPriorDialog implements AbstractPriorDialog {

    private static final int MINIMUM_TABLE_WIDTH = 120;
    private static final int MINIMUM_TABLE_HEIGHT = 160;
    private static final int PREFERRED_TABLE_WIDTH = 180;
    private static final int PREFERRED_TABLE_HEIGHT = 320;

    private JFrame frame;

    private final JTable parametersTable;
    private final ParametersTableModel parametersTableModel;

    private JTextField nameField = new JTextField();
    private JPanel panel;

    private final PriorSettingsPanel priorSettingsPanel;

    private Parameter parameter;

    final private BeautiOptions options;

    private SelectParametersDialog selectParametersDialog = null;
    private List<Parameter> compatibleParameterList;
    private List<Parameter> dependentParameterList;

    public JointPriorDialog(JFrame frame, BeautiOptions options) {
        this.frame = frame;

        this.options = options;

        priorSettingsPanel = new PriorSettingsPanel(frame);

        nameField.setColumns(30);

        parametersTableModel = new ParametersTableModel();
//        TableSorter sorter = new TableSorter(traitsTableModel);
//        traitsTable = new JTable(sorter);
//        sorter.setTableHeader(traitsTable.getTableHeader());
        parametersTable = new JTable(parametersTableModel);

        parametersTable.getTableHeader().setReorderingAllowed(false);
        parametersTable.getTableHeader().setResizingAllowed(false);
//        traitsTable.getTableHeader().setDefaultRenderer(
//                new HeaderRenderer(SwingConstants.LEFT, new Insets(0, 4, 0, 4)));

        parametersTable.getSelectionModel().setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        parametersTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent evt) {
                parametersSelectionChanged();
            }
        });

    }

    private void parametersSelectionChanged() {
        int selRow = parametersTable.getSelectedRow();
        if (selRow >= 0) {
            removeParameterAction.setEnabled(true);
        }

        if (dependentParameterList.size() <= 0) {
            removeParameterAction.setEnabled(false);
        }
    }

    public boolean validateModelName() {
        return validateModelName(nameField.getText());
    }

    private boolean validateModelName(String name) {
//        System.err.println("Validating: " + modelName);
        // check that the name is valid

        if (name.equals(parameter.getName())) {
            return true;
        }

        if (name.trim().length() == 0) {
            Toolkit.getDefaultToolkit().beep();
            return false;
        }

        // check that a parameter with this name doesn't exist
        if (parameterExists(name)) {
            JOptionPane.showMessageDialog(frame,
                    "A parameter with this name already exists.",
                    "Linked parameter error",
                    JOptionPane.WARNING_MESSAGE);
            return false;
        }

        // check that a model with this name doesn't exist
        if (modelExists(name)) {
            JOptionPane.showMessageDialog(frame,
                    "A model with this name already exists.",
                    "Linked parameter error",
                    JOptionPane.WARNING_MESSAGE);
            return false;
        }

        return true;
    }

    private boolean parameterExists(String name) {
        for (Parameter parameter : options.selectParameters()) {
            if (parameter.getName().equals(name)) {
                return true;
            }
        }
        return false;
    }

    private boolean modelExists(String modelName) {
        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);
        return comp.modelExists(modelName);
    }

    public int showDialog() {
        panel = new JPanel(new GridBagLayout());

        double lower = Double.NEGATIVE_INFINITY;
        double upper = Double.POSITIVE_INFINITY;

        if (parameter.isZeroOne) {
            lower = 0.0;
            upper = 1.0;
        } else if (parameter.isNonNegative) {
            lower = 0.0;
        }

        panel = new JPanel(new GridBagLayout());

        setupComponents();

        JOptionPane optionPane = new JOptionPane(panel,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Linked Parameter Setup");

        priorSettingsPanel.setDialog(dialog);
        priorSettingsPanel.setParameter(parameter);

        if (OSType.isMac()) {
            dialog.setMinimumSize(new Dimension(dialog.getBounds().width, 300));
        } else {
            Toolkit tk = Toolkit.getDefaultToolkit();
            Dimension d = tk.getScreenSize();
            if (d.height < 700 && panel.getHeight() > 450) {
                dialog.setSize(new Dimension(panel.getWidth() + 100, 550));
            } else {
                // setSize because optionsPanel is shrunk in dialog
                dialog.setSize(new Dimension(panel.getWidth() + 100, panel.getHeight() + 100));
            }

//            System.out.println("panel width = " + panel.getWidth());
//            System.out.println("panel height = " + panel.getHeight());
        }

        dialog.pack();
        dialog.setResizable(true);
        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public String getName() {
        return nameField.getText();
    }

    public List<Parameter> getDependentParameterList() {
        return dependentParameterList;
    }

    private void setupComponents() {
        panel.removeAll();

        JScrollPane scrollPane1 = new JScrollPane(parametersTable,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane1.setOpaque(false);

        ActionPanel actionPanel1 = new ActionPanel(false);
        actionPanel1.setAddAction(addParameterAction);
        actionPanel1.setRemoveAction(removeParameterAction);
        actionPanel1.setAddToolTipText("Use this button to add an existing parameter to the prior");
        actionPanel1.setRemoveToolTipText("Use this button to remove a parameter from the prior");

        removeParameterAction.setEnabled(false);

        JPanel controlPanel1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
        controlPanel1.setOpaque(false);
        controlPanel1.add(actionPanel1);

        JPanel panel1 = new JPanel(new BorderLayout(0, 0));
        panel1.setOpaque(false);
        panel1.add(new JLabel("Linked parameters:"), BorderLayout.NORTH);
        panel1.add(scrollPane1, BorderLayout.CENTER);

        // removing the control panel for now. Not sure whether we really want adding and
        // removing of parameteres in this dialog.
//        panel1.add(controlPanel1, BorderLayout.SOUTH);
        panel1.setSize(new Dimension(PREFERRED_TABLE_WIDTH, PREFERRED_TABLE_HEIGHT));
        panel1.setPreferredSize(new Dimension(PREFERRED_TABLE_WIDTH, PREFERRED_TABLE_HEIGHT));
        panel1.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, MINIMUM_TABLE_HEIGHT));

        OptionsPanel optionsPanel = new OptionsPanel(0,6);
        if (parameter.getName() != null) {
            nameField.setText(parameter.getName());
        } else {
            nameField.setText("Untitled");
        }
        optionsPanel.addComponentWithLabel("Unique Name: ", nameField);
//        optionsPanel.addComponentWithLabel("Initial Value: ", initialField);

        panel.setOpaque(false);
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));

        GridBagConstraints c = new GridBagConstraints();
        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0;
        c.weighty = 0;
        c.anchor = GridBagConstraints.PAGE_START;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridwidth = 2;
        panel.add(optionsPanel, c);

        c.gridx = 0;
        c.gridy = 1;
        c.anchor = GridBagConstraints.PAGE_START;
        c.fill = GridBagConstraints.VERTICAL;
        c.gridwidth = 1;
        panel.add(panel1, c);

        c.gridx = 1;
        c.gridy = 1;
        c.weightx = 1;
        c.weighty = 1;
        c.fill = GridBagConstraints.BOTH;
        c.anchor = GridBagConstraints.PAGE_START;
        c.gridwidth = GridBagConstraints.REMAINDER;
        panel.add(priorSettingsPanel, c);

    }

    public void getArguments(Parameter parameter) {
        priorSettingsPanel.getArguments(parameter);
    }

    public boolean hasInvalidInput(boolean showError) {
        return priorSettingsPanel.hasInvalidInput(showError);
    }

    public boolean addParameter() {
        if (selectParametersDialog == null) {
            selectParametersDialog = new SelectParametersDialog(frame);
        }
        List<Parameter> availableParameters = new ArrayList<Parameter>(compatibleParameterList);
        availableParameters.removeAll(dependentParameterList);
        int result = selectParametersDialog.showDialog("Select parameter to add to this Linked Parameter", availableParameters);
        if (result == JOptionPane.OK_OPTION) {
            Parameter parameter = selectParametersDialog.getSelectedParameter();
            dependentParameterList.add(parameter);
            parametersTableModel.fireTableDataChanged();
        } else if (result == JOptionPane.CANCEL_OPTION) {
            return false;
        }

        return true;
    }

    private void removeSelectedParameters() {
        int[] selRows = parametersTable.getSelectedRows();
        List<Parameter> parametersToRemove = new ArrayList<Parameter>();
        for (int row : selRows) {
            parametersToRemove.add((Parameter)parametersTable.getValueAt(row, 0));
        }
        removeParameters(parametersToRemove);
    }

    private void removeParameters(List<Parameter> parametersToRemove) {
        for (Parameter parameter : parametersToRemove) {
            dependentParameterList.remove(parameter);
        }
        parametersTableModel.fireTableDataChanged();
    }

    private AddParameterAction addParameterAction = new AddParameterAction();

    public void setLinkedParameter(LinkedParameter linkedParameter) {
        parameter = linkedParameter.getArgumentParameter();
    }

    public void setDependentParameterList(List<Parameter> dependentParameterList) {
        this.dependentParameterList = dependentParameterList;
    }

    public void setCompatibleParameterList(List<Parameter> compatibleParameterList) {
        this.compatibleParameterList = compatibleParameterList;
    }

    public class AddParameterAction extends AbstractAction {

        public AddParameterAction() {
            super("Add parameter");
        }

        public void actionPerformed(ActionEvent ae) {
            addParameter();
        }
    }

    AbstractAction removeParameterAction = new AbstractAction() {
        public void actionPerformed(ActionEvent ae) {
            removeSelectedParameters();
        }
    };


    class ParametersTableModel extends AbstractTableModel {

        private static final long serialVersionUID = -6707994233020715574L;
        String[] columnNames = {"Parameter"};

        public ParametersTableModel() {
        }

        public int getColumnCount() {
            return columnNames.length;
        }

        public int getRowCount() {
            if (dependentParameterList == null) {
                return 0;
            }
            return dependentParameterList.size();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return dependentParameterList.get(row);
            }
            return null;
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