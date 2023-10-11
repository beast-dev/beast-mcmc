package dr.app.beauti.datapanel;

import dr.app.beauti.BeautiFrame;
import dr.app.beauti.traitspanel.CreateTraitDialog;
import jam.panels.OptionsPanel;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class CreateBadTraitFormatDialog {

    private final BeautiFrame frame;
    private final JButton exampleButton = new JButton("Show example of mapping file format");



    OptionsPanel optionPanel;
    private final JLabel label = new JLabel();

    public CreateBadTraitFormatDialog(final BeautiFrame frame) {
        this.frame = frame;
        optionPanel = new OptionsPanel(12, 12);

        optionPanel.addSpanningComponent(new JLabel("An error occurred when importing the traits file."));
        optionPanel.addSpanningComponent(new JLabel("It is likely the file is not properly formatted."));

        optionPanel.addComponent(exampleButton);

        exampleButton.setEnabled(true);
        exampleButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ae) {
                CreateTraitDialog.showExampleTraitFormat(frame);
            }
        });
    }

    public int showDialog() {
        String[] options = {"OK"};

        int choice = JOptionPane.showOptionDialog(frame,
                optionPanel,
                "Could not import traits.",
                JOptionPane.DEFAULT_OPTION,
                JOptionPane.ERROR_MESSAGE,
                null,
                options,
                options[0]);
        return choice;
    }

}
