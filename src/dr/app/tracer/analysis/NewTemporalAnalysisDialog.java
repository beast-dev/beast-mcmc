package dr.app.tracer.analysis;

import org.virion.jam.components.RealNumberField;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.framework.DocumentFrame;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;

public class NewTemporalAnalysisDialog {

	private JFrame frame;

	private JTextField titleField;

	private WholeNumberField binCountField;

	private RealNumberField minTimeField;
	private RealNumberField maxTimeField;

	private OptionsPanel optionPanel;

	public NewTemporalAnalysisDialog(JFrame frame) {
		this.frame = frame;

		titleField = new JTextField();
		titleField.setColumns(32);

		binCountField = new WholeNumberField(2, 2000);
		binCountField.setValue(100);
		binCountField.setColumns(4);

		maxTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
		maxTimeField.setColumns(12);

		minTimeField = new RealNumberField(0.0, Double.MAX_VALUE);
		minTimeField.setColumns(12);

		optionPanel = new OptionsPanel(12, 12);
	}

	public int showDialog() {

		setArguments();

		JOptionPane optionPane = new JOptionPane(optionPanel,
				JOptionPane.QUESTION_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				null,
				null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

		final JDialog dialog = optionPane.createDialog(frame, "Demographic Analysis");
		dialog.pack();

		dialog.setVisible(true);

		int result = JOptionPane.CANCEL_OPTION;
		Integer value = (Integer)optionPane.getValue();
		if (value != null && value.intValue() != -1) {
			result = value.intValue();
		}

		if (result == JOptionPane.OK_OPTION) {
		}

		return result;
	}

	private void setArguments() {
		optionPanel.removeAll();

		optionPanel.addComponentWithLabel("Title:", titleField);

		optionPanel.addSeparator();

		optionPanel.addComponentWithLabel("Number of bins:", binCountField);

		optionPanel.addComponentWithLabel("Minimum time:", minTimeField);
		optionPanel.addComponentWithLabel("Maximum time:", maxTimeField);

//		JLabel label3 = new JLabel(
//				"<html>You can set the age of sampling of the most recent tip in<br>" +
//						"the tree. If this is set to zero then the plot is shown going<br>" +
//						"backwards in time, otherwise forwards in time.</html>");
//		label3.setFont(label3.getFont().deriveFont(((float)label3.getFont().getSize() - 2)));
//		optionPanel.addSpanningComponent(label3);
	}

	public TemporalAnalysisFrame createTemporalAnalysisFrame(DocumentFrame parent) {

		TemporalAnalysisFrame frame = new TemporalAnalysisFrame(parent, titleField.getText(),
				binCountField.getValue(), minTimeField.getValue(), maxTimeField.getValue());
		frame.initialize();
		return frame;
	}

}
