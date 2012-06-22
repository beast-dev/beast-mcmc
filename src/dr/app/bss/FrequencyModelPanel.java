package dr.app.bss;

import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingConstants;

import org.virion.jam.components.RealNumberField;

@SuppressWarnings("serial")
public class FrequencyModelPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JScrollPane scrollPane;
	private JComboBox frequencyCombo;
	private RealNumberField[] frequencyParameterFields = new RealNumberField[BeagleSequenceSimulatorData.frequencyParameterNames.length];
	
	public FrequencyModelPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {
		
		this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		scrollPane = new JScrollPane();
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		scrollPane = new JScrollPane(optionPanel,
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollPane.setOpaque(false);
		
		add(scrollPane, "Center");
		
		frequencyCombo = new JComboBox();
		frequencyCombo.setOpaque(false);

		for (String frequencyModel : BeagleSequenceSimulatorData.frequencyModels) {
			frequencyCombo.addItem(frequencyModel);
		}// END: fill loop

		frequencyCombo.addItemListener(new ListenFrequencyCombo());

		for (int i = 0; i < BeagleSequenceSimulatorData.frequencyParameterNames.length; i++) {
			frequencyParameterFields[i] = new RealNumberField();
			frequencyParameterFields[i].setColumns(8);
			frequencyParameterFields[i].setValue(data.frequencyParameterValues[i]);
		}// END: fill loop

		setFrequencyArguments();
		
	}//END: Constructor

	private void setFrequencyArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Frequency model:"), frequencyCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = frequencyCombo.getSelectedIndex();

		for (int i = 0; i < data.frequencyParameterIndices[index].length; i++) {

			int k = data.frequencyParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(frequencyParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					BeagleSequenceSimulatorData.frequencyParameterNames[k] + ":",
					panel);

		}// END: indices loop

//		validate();
//		repaint();
	}// END: setFrequencyArguments

	private class ListenFrequencyCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setFrequencyArguments();
			frame.fireModelChanged();

		}// END: actionPerformed
	}// END: ListenClockCombo

	public void collectSettings() {

		data.frequencyModel = frequencyCombo.getSelectedIndex();
		for (int i = 0; i < BeagleSequenceSimulatorData.frequencyParameterNames.length; i++) {

			data.frequencyParameterValues[i] = frequencyParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}//END: getExportableComponent
	
}//END: class
