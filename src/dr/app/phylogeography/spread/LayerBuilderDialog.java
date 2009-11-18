package dr.app.phylogeography.spread;

import dr.app.phylogeography.builder.BuilderFactory;
import dr.app.phylogeography.builder.Builder;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class LayerBuilderDialog {

    private JFrame frame;

    private JComboBox builderCombo;
    private JTextField nameField;

    private OptionsPanel optionPanel;

    private Builder builder = null;

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
        if (builder != null) {
            optionPanel.addSeparator();
            optionPanel.addSpanningComponent(builder.getEditPanel());
        }
    }

    public int showDialog(BuilderFactory[] builderFactories) {

        builderCombo = new JComboBox();
        builderCombo.addItem("Select Layer type...");
        for (Object factory : builderFactories) {
            builderCombo.addItem(factory);
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

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public int showDialog(Builder builder) {

        builderCombo = null;
        
        this.builder = builder;
        nameField.setText(builder.getName());

        setupPanel(builder);

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Edit Layer Settings");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public Builder getBuilder() {
        builder.setName(nameField.getText());
        builder.setFromEditPanel();
        return builder;
    }

}