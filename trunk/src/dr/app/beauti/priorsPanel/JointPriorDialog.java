/*
 * PriorDialog.java
 *
 * @author Marc A. Suchard
 */

package dr.app.beauti.priorsPanel;

import dr.app.beauti.ComboBoxRenderer;
import dr.app.beauti.components.hpm.HierarchicalModelComponentOptions;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.traitspanel.CreateTraitDialog;
import dr.app.beauti.types.PriorType;
import dr.app.gui.components.RealNumberField;
import dr.app.gui.table.TableEditorStopper;
import dr.app.util.OSType;
import jam.panels.ActionPanel;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.BorderUIResource;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableColumn;
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
public class JointPriorDialog {

    private static final int MINIMUM_TABLE_WIDTH = 140;

    private JFrame frame;

    private final JTable parametersTable;
    private final ParametersTableModel parametersTableModel;

    private RealNumberField initialField = new RealNumberField();
    private JTextField nameField = new JTextField();
    private JPanel panel;

    private java.util.List<Parameter> parameterList;
    private Parameter parameter;

    final private BeautiOptions options;

    public JointPriorDialog(JFrame frame, BeautiOptions options) {
        this.frame = frame;
        this.options = options;

        nameField.setColumns(10);
        initialField.setColumns(10);

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

        if (parameterList.size() <= 0) {
            removeParameterAction.setEnabled(false);
        }
    }


    public void addParameters(java.util.List<Parameter> parameterList) {
    }

    public boolean validateModelName() {
        return validateModelName(nameField.getText());
    }

    private boolean validateModelName(String modelName) {
//        System.err.println("Validating: " + modelName);
        // check that the name is valid
        if (modelName.trim().length() == 0) {
            Toolkit.getDefaultToolkit().beep();
            return false;
        }

        // check that the trait name doesn't exist
        if (modelExists(modelName)) {
            JOptionPane.showMessageDialog(frame,
                    "A model with this name already exists.",
                    "HPM name error",
//                    JOptionPane.YES_NO_OPTION,
                    JOptionPane.WARNING_MESSAGE);
//            System.err.println("Model name exists");
            return false;

//            if (option == JOptionPane.NO_OPTION) {
//                return false;
//            }
        }

        return true;
    }

    private boolean modelExists(String modelName) {

        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);
        return comp.modelExists(modelName);
    }

    public int showDialog(final java.util.List<Parameter> parameterList) {

        this.parameterList = parameterList;
        this.parameter = parameterList.get(0);

        panel = new JPanel(new GridBagLayout());

        double lower = Double.NEGATIVE_INFINITY;
        double upper = Double.POSITIVE_INFINITY;

        if (parameter.isZeroOne) {
            lower = 0.0;
            upper = 1.0;
        } else if (parameter.isNonNegative) {
            lower = 0.0;
        }

        initialField.setRange(lower, upper);
        initialField.setValue(parameter.initial);

        panel = new JPanel(new GridBagLayout());

        setupComponents();

        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        JOptionPane optionPane = new JOptionPane(scrollPane,
                JOptionPane.PLAIN_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Phylogenetic Hierarchical Model Setup");


        dialog.pack();
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

        dialog.setResizable(true);
        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

//        if (result == JOptionPane.OK_OPTION) {
//            getArguments();
//        }

        return result;
    }

    public String getName() {
        return nameField.getText();
    }

    public List<Parameter> getParameterList() {
        return parameterList;
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
        panel1.add(scrollPane1, BorderLayout.CENTER);
        panel1.add(controlPanel1, BorderLayout.SOUTH);
        panel1.setMinimumSize(new Dimension(MINIMUM_TABLE_WIDTH, 0));

        OptionsPanel optionsPanel = new OptionsPanel(0,6);
        optionsPanel.addComponentWithLabel("Unique Name: ", nameField);
        optionsPanel.addComponentWithLabel("Initial Value: ", initialField);

        //panel2.add(toolBar1, BorderLayout.NORTH);
        //panel2.add(scrollPane2, BorderLayout.CENTER);

//        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
//                panel1, panel2);
//        splitPane.setDividerLocation(MINIMUM_TABLE_WIDTH);
//        splitPane.setContinuousLayout(true);
//        splitPane.setBorder(BorderFactory.createEmptyBorder());
//        splitPane.setOpaque(false);

        panel.setOpaque(false);
        panel.setBorder(new BorderUIResource.EmptyBorderUIResource(new Insets(12, 12, 12, 12)));
        panel.setLayout(new BorderLayout(0, 0));
        panel.add(optionsPanel, BorderLayout.CENTER);
        panel.add(panel1, BorderLayout.EAST);
        //panel.add(toolBar1, BorderLayout.NORTH);
    }

    protected JPanel createPriorPanel() {
        return new JPanel();
    }

    public boolean addParameter() {
//        if (addParametersDialog == null) {
//            addParametersDialog = new AddParametersDialog(frame);
//        }
//
//        addParametersDialog.setSpeciesTrait(isSpeciesTrait);
//        addParametersDialog.setTraitName(traitName);
//        addParametersDialog.setMessage(message);
//
//        int result = addParametersDialog.showDialog();
//        if (result == JOptionPane.OK_OPTION) {
//            fireParametersChanged();
//            updateButtons();
//
//        } else if (result == JOptionPane.CANCEL_OPTION) {
//            return false;
//        }

        return true;
    }

    private void removeParameter() {
        int selRow = parametersTable.getSelectedRow();
        removeParameter((Parameter)parametersTable.getValueAt(selRow, 0));
    }

    private void removeParameter(Parameter parameter) {
    }

    private AddParameterAction addParameterAction = new AddParameterAction();

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
            removeParameter();
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
            if (parameterList == null) {
                return 0;
            }
            return parameterList.size();
        }

        public Object getValueAt(int row, int col) {
            switch (col) {
                case 0:
                    return parameterList.get(row).getName();
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