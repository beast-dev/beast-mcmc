package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.virion.jam.components.WholeNumberField;

@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

	// private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JLabel replicatesLabel = new JLabel("Number of replicates:");
	private WholeNumberField replicatesField = new WholeNumberField(1,
			Integer.MAX_VALUE);

	public SimulationPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {

		// this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		replicatesField.setColumns(8);
		replicatesField.setValue(data.replicateCount);

		optionPanel.addComponents(replicatesLabel, replicatesField);

	}// END: SimulationPanel

	public final void collectSettings() {
		data.replicateCount = replicatesField.getValue();
	}

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
