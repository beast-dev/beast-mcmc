package pyromania.sourcespanel;

import jam.panels.OptionsPanel;
import pyromania.database.Source;
import pyromania.platforms.PlatformType;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class SourcesSettingsDialog {

    private JFrame frame;

    private OptionsPanel optionPanel;

    private Source source = null;

    private JTextField nameField = new JTextField();
    private JTextField descriptionField = new JTextField();
    private JComboBox platformCombo = new JComboBox(PlatformType.values());

    public SourcesSettingsDialog(JFrame frame) {
        this.frame = frame;

        nameField.setColumns(32);
        descriptionField.setColumns(32);

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel(Source source) {
        optionPanel.removeAll();

        nameField.setText(source.getName());
        optionPanel.addComponentWithLabel("Name:", nameField);

        descriptionField.setText(source.getDescription());
        optionPanel.addComponentWithLabel("Description:", descriptionField);

        platformCombo.setSelectedItem(source.getPlatformType());
        optionPanel.addComponentWithLabel("Platform:", platformCombo);

    }

    public int showDialog(Source source) {

        this.source = source;

        setupPanel(source);

        final JOptionPane optionPane = new JOptionPane(optionPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, "Edit Reads File Settings");
        dialog.pack();


        int result = JOptionPane.CANCEL_OPTION;

        dialog.setVisible(true);

        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

    public Source getReadsFile() {
        source.setName(nameField.getText());
        source.setDescription(descriptionField.getText());
        source.setPlatformType((PlatformType)platformCombo.getSelectedItem());
        return source;
    }

}