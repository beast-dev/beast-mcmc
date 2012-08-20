package dr.app.tracer.traces;

import javax.swing.*;
import java.awt.*;

/**
 * @author Walter Xie
 */
public class FilterContinuousPanel extends FilterAbstractPanel {
    JTextField minField;
    JTextField maxField;

    FilterContinuousPanel(String[] minMax, String[] bound) {
        setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();

        if (bound == null) {
            bound = new String[2];
        }

        minField = new JTextField(bound[0]);
        minField.setColumns(20);

//        c.weightx = 5;
//        c.weighty = 10;
//        c.gridx = 0;
        c.gridy = 0;
        c.insets = new Insets(20,10,0,10);
        c.anchor = GridBagConstraints.FIRST_LINE_START;
        add(new JLabel("Set Minimum for Selecting Values : "), c);

        c.gridy = 1;
        add(minField, c);

        c.gridy = 2;
        add(new JLabel("which should > " + minMax[0]), c);

        maxField = new JTextField(bound[1]);
        maxField.setColumns(20);

        c.gridy = 3;
        c.insets = new Insets(50,10,0,10);
        add(new JLabel("Set Maximum for Selecting Values : "), c);

        c.gridy = 4;
        c.insets = new Insets(20,10,0,10);
        add(maxField, c);

        c.gridy = 5;   
        add(new JLabel("which should < " + minMax[1]), c);
    }

    public Object[] getSelectedValues() {
        return new String[]{minField.getText(), maxField.getText()};
    }

}
