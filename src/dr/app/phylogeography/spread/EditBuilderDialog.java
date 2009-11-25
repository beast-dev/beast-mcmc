package dr.app.phylogeography.spread;

import dr.app.phylogeography.builder.Builder;
import dr.app.phylogeography.builder.BuilderFactory;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

/**
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class EditBuilderDialog {

    private JFrame frame;

    private JTextField nameField;

    private OptionsPanel optionPanel;

    private Builder builder = null;

    public EditBuilderDialog(JFrame frame) {
        this.frame = frame;

        nameField = new JTextField();
        nameField.setColumns(20);

        optionPanel = new OptionsPanel(12, 12);
    }

    private void setupPanel(Builder builder) {
        optionPanel.removeAll();

        optionPanel.addComponentWithLabel("Name:", nameField);

        optionPanel.addSeparator();
        optionPanel.addSpanningComponent(builder.getEditPanel());
    }

    public int showDialog(Builder builder, SpreadDocument document) {

        this.builder = builder;
        nameField.setText(builder.getName());

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
                if (nameField.getText().trim().length() > 0) {
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
        builder.setFromEditPanel();
        return builder;
    }

}