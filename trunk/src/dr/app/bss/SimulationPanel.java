package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.virion.jam.components.WholeNumberField;

@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

//	 private BeagleSequenceSimulatorFrame frame;
	private PartitionDataList dataList;
	private OptionsPanel optionPanel;

	private JLabel replicatesLabel = new JLabel("Number of sites:");
	private WholeNumberField replicatesField = new WholeNumberField(1,
			Integer.MAX_VALUE);
	
	// Buttons
	private JButton simulate;
	
	public SimulationPanel(
//			final BeagleSequenceSimulatorFrame frame,
			final PartitionDataList dataList) {

//		this.frame = frame;
		this.dataList = dataList;

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		
		// number of sites
		replicatesField.setColumns(8);
		replicatesField.setValue(dataList.sequenceLength);
		
		optionPanel.addComponents(replicatesLabel, replicatesField);

		// Buttons
		simulate = new JButton("Simulate", BeagleSequenceSimulatorApp.hammerIcon);
		simulate.addActionListener(new ListenSimulate());
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);
		buttonsHolder.add(simulate);

		setOpaque(false);
		setLayout(new BorderLayout());
		add(optionPanel, BorderLayout.NORTH);
		add(buttonsHolder, BorderLayout.SOUTH);

	}// END: SimulationPanel

	public final void collectSettings() {

		// all elements hold the same value
			dataList.sequenceLength = replicatesField.getValue();
		// frame.fireModelChanged();
	}// END: collectSettings

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	private class ListenSimulate implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

//			frame.doExport();
			Utils.printDataList(dataList);
			
		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates
	
}// END: class
