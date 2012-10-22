package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.util.ArrayList;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.virion.jam.components.WholeNumberField;

@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

	 private BeagleSequenceSimulatorFrame frame;
	private ArrayList<PartitionData> dataList;
	private OptionsPanel optionPanel;

	private JLabel replicatesLabel = new JLabel("Number of sites:");
	private WholeNumberField replicatesField = new WholeNumberField(1,
			Integer.MAX_VALUE);
	
	private JPanel panel;
	
	public SimulationPanel(final BeagleSequenceSimulatorFrame frame,
			final ArrayList<PartitionData> dataList) {

		this.frame = frame;
		this.dataList = dataList;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		// number of sites
		replicatesField.setColumns(8);
		replicatesField.setValue(dataList.get(0).sequenceLength);
		
		optionPanel.addComponents(replicatesLabel, replicatesField);

		// simulation type
		panel = new JPanel();
		panel.setLayout(new GridLayout(2, 1));
		
	}// END: SimulationPanel

	public final void collectSettings() {
		dataList.get(0).sequenceLength = replicatesField.getValue();
		frame.fireModelChanged();
	}

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
