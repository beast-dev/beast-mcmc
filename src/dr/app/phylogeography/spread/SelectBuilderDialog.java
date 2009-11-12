package dr.app.phylogeography.spread;

import dr.app.phylogeography.builder.BuilderFactory;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class SelectBuilderDialog {

    private JFrame frame;

    JComboBox builderCombo;
    JTextField nameField;

    OptionsPanel optionPanel;

    public SelectBuilderDialog(JFrame frame) {
        this.frame = frame;

        builderCombo = new JComboBox();

        nameField = new JTextField();
        nameField.setColumns(20);

        optionPanel = new OptionsPanel(12, 12);
        optionPanel.addComponentWithLabel("Name:", nameField);
        optionPanel.addComponentWithLabel("Layer type:", builderCombo);

    }

    public int showDialog(BuilderFactory[] builderFactories) {

        builderCombo.removeAllItems();
        for (Object factory : builderFactories) {
            builderCombo.addItem(factory);
        }

        JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Create New Partition Clock Model");
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public BuilderFactory getBuilderFactory() {
        return (BuilderFactory) builderCombo.getSelectedItem();
    }

    public String getName() {
        return nameField.getText();
    }

}