package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import dr.app.gui.components.RealNumberField;

@SuppressWarnings("serial")
public class ModelsPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JComboBox substitutionCombo;
	private RealNumberField[] substitutionParameterFields = new RealNumberField[BeagleSequenceSimulatorData.substitutionParameterNames.length];

	public ModelsPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {

		super();

		this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);
		
		substitutionCombo = new JComboBox();
		substitutionCombo.setOpaque(false);

		for (String substitutionModel : BeagleSequenceSimulatorData.substitutionModels) {
			substitutionCombo.addItem(substitutionModel);
		}// END: fill loop

		substitutionCombo.addItemListener(new ListenSubstitutionCombo());

		for (int i = 0; i < BeagleSequenceSimulatorData.substitutionParameterNames.length; i++) {
			substitutionParameterFields[i] = new RealNumberField();
			substitutionParameterFields[i].setColumns(8);
			substitutionParameterFields[i].setValue(data.substitutionParameterValues[i]);
		}// END: fill loop

		setSubstitutionArguments();
		
		
	}// END: Constructor

	private void setSubstitutionArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Branch substitution model:"), substitutionCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Select substitution parameter values:");

		int substIndex = substitutionCombo.getSelectedIndex();

		for (int i = 0; i < data.substitutionParameterIndices[substIndex].length; i++) {

			int k = data.substitutionParameterIndices[substIndex][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(substitutionParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(BeagleSequenceSimulatorData.substitutionParameterNames[k] + ":", panel);

		}// END: indices loop

		validate();
		repaint();
	}// END: setSubstitutionArguments

	private class ListenSubstitutionCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setSubstitutionArguments();
			frame.fireModelChanged();

		}// END: actionPerformed
	}// END: ListenSubstitutionCombo

	public void collectSettings() {

		data.substitutionModel = substitutionCombo.getSelectedIndex();
		for (int i = 0; i < BeagleSequenceSimulatorData.substitutionParameterNames.length; i++) {

			data.substitutionParameterValues[i] = substitutionParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	public JComponent getExportableComponent() {
		return this;
	}//END: getExportableComponent

}// END: class
