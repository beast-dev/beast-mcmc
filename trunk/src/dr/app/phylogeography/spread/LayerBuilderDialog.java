package dr.app.phylogeography.spread;

import dr.app.phylogeography.builder.BuilderFactory;
import dr.app.phylogeography.builder.Builder;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.text.Document;
import javax.swing.border.EmptyBorder;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class LayerBuilderDialog {

    private JFrame frame;

    private JComboBox dataFileCombo;
    private JComboBox builderCombo;
    private JTextField nameField;

    private OptionsPanel optionPanel;

    private Builder builder = null;

    private static final String SELECT_LAYER_TYPE = "Select layer type...";
    private static final String SELECT_INPUT_FILE = "Select input file...";

    public LayerBuilderDialog(JFrame frame) {
        this.frame = frame;

        nameField = new JTextField();
        nameField.setColumns(20);

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel(Builder builder) {
        optionPanel.removeAll();

        optionPanel.addComponentWithLabel("Name:", nameField);

        if (builderCombo != null) {
            optionPanel.addComponentWithLabel("Layer type:", builderCombo);
        }

        if (dataFileCombo != null) {
            optionPanel.addComponentWithLabel("Input file:", dataFileCombo);
        }

        if (builder != null) {
            optionPanel.addSeparator();
            optionPanel.addSpanningComponent(builder.getEditPanel());
        }
    }

    public int showDialog(BuilderFactory[] builderFactories, SpreadDocument document) {

        builderCombo = new JComboBox();
        builderCombo.addItem(SELECT_LAYER_TYPE);
        for (Object factory : builderFactories) {
            builderCombo.addItem(factory);
        }

        dataFileCombo = new JComboBox();
        dataFileCombo.addItem(SELECT_INPUT_FILE);
        for (Object dataFile : document.getDataFiles()) {
            dataFileCombo.addItem(dataFile);
        }

        setupPanel(null);

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Layer");
        dialog.pack();

        builderCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent e) {
                if (builderCombo.getSelectedIndex() != 0) {
                    builder = ((BuilderFactory)builderCombo.getSelectedItem()).createBuilder();
                    setupPanel(builder);
                } else {
                    setupPanel(null);
                }
                dialog.pack();
            }
        });

        return validateDialog(dialog, optionPane);
    }

    public int showDialog(Builder builder, SpreadDocument document) {

        builderCombo = null;

        this.builder = builder;
        nameField.setText(builder.getName());

        dataFileCombo = new JComboBox();
        for (Object dataFile : document.getDataFiles()) {
            dataFileCombo.addItem(dataFile);
        }
        dataFileCombo.setSelectedItem(builder.getDataFile());

        setupPanel(builder);

        final JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Edit Layer Settings");
        dialog.pack();

        return validateDialog(dialog, optionPane);
    }

    private int validateDialog(final JDialog dialog, final JOptionPane optionPane) {

        int result;
        boolean validated = false;

        do {
            result = JOptionPane.CANCEL_OPTION;

            dialog.setVisible(true);

            Integer value = (Integer) optionPane.getValue();
            if (value != null && value != -1) {
                result = value;
            }

            if (result != JOptionPane.CANCEL_OPTION) {
                if (nameField.getText().trim().length() > 0 &&
                        (dataFileCombo == null || !dataFileCombo.getSelectedItem().toString().equals(SELECT_LAYER_TYPE)) &&
                        (builderCombo == null || !builderCombo.getSelectedItem().toString().equals(SELECT_INPUT_FILE))) {
                    validated = true;
                } else {
                    Toolkit.getDefaultToolkit().beep();
                }
            }

        } while (!validated && result != JOptionPane.CANCEL_OPTION);

        return result;
    }

    public Builder getBuilder() {
        builder.setName(nameField.getText());
        builder.setDataFile((SpreadDocument.DataFile)dataFileCombo.getSelectedItem());
        builder.setFromEditPanel();
        return builder;
    }

}