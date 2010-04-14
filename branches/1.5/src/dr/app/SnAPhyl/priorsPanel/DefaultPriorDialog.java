package dr.app.SnAPhyl.priorsPanel;

import dr.app.SnAPhyl.SnAPhylFrame;
import dr.app.beauti.options.BeautiOptions;

import javax.swing.*;
import javax.swing.border.EmptyBorder;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: PriorDialog.java,v 1.4 2006/09/05 13:29:34 rambaut Exp $
 */
public class DefaultPriorDialog {

   private SnAPhylFrame frame;
   private PriorsPanel priorsPanel;
    
	public DefaultPriorDialog(SnAPhylFrame frame) {
		this.frame = frame;
        priorsPanel = new PriorsPanel(frame, true);
	}

	public boolean showDialog(BeautiOptions options) {
        priorsPanel.setParametersList(options);

        Object[] buttons = {"Continue", "Back to BEAUTi"};
        JOptionPane optionPane = new JOptionPane(priorsPanel,
				JOptionPane.PLAIN_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				buttons,
				buttons[0]);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));
//        optionPane.setPreferredSize(new java.awt.Dimension(800, 600));

		final JDialog dialog = optionPane.createDialog(frame, "Check the default priors");
		dialog.pack();
		dialog.setVisible(true);

        return optionPane.getValue() != null && optionPane.getValue().equals(buttons[0]);

    }
}