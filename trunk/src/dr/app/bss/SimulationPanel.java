package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
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

	private BeagleSequenceSimulatorFrame frame;
	private PartitionDataList dataList;
	private OptionsPanel optionPanel;

	private JLabel sequenceLengthLabel;
	private WholeNumberField sequenceLengthField;

	private JLabel simulationsNumberLabel;
	private WholeNumberField simulationsNumberField;

	// Buttons
	private JButton simulate;
	private JButton generateXML;

	public SimulationPanel(final BeagleSequenceSimulatorFrame frame,
			final PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		// number of simulations
		simulationsNumberLabel = new JLabel("Number of simulations:");
		simulationsNumberField = new WholeNumberField(1, Integer.MAX_VALUE);
		simulationsNumberField.setColumns(8);
		simulationsNumberField.setValue(dataList.simulationsNumber);
		optionPanel.addComponents(simulationsNumberLabel,
				simulationsNumberField);

		// number of sites
		sequenceLengthLabel = new JLabel("Number of sites:");
		sequenceLengthField = new WholeNumberField(1, Integer.MAX_VALUE);
		sequenceLengthField.setColumns(8);
		sequenceLengthField.setValue(dataList.sequenceLength);
		optionPanel.addComponents(sequenceLengthLabel, sequenceLengthField);

		// Buttons holder
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);

		// simulate button
		simulate = new JButton("Simulate",
				BeagleSequenceSimulatorApp.nuclearIcon);
		simulate.addActionListener(new ListenSimulate());
		buttonsHolder.add(simulate);

		generateXML = new JButton("Generate XML",
				BeagleSequenceSimulatorApp.hammerIcon);
		generateXML.addActionListener(new ListenGenerateXML());
		buttonsHolder.add(generateXML);

		setOpaque(false);
		setLayout(new BorderLayout());
		add(optionPanel, BorderLayout.NORTH);
		add(buttonsHolder, BorderLayout.SOUTH);

	}// END: SimulationPanel

	public final void collectSettings() {
		dataList.sequenceLength = sequenceLengthField.getValue();
		// frame.fireModelChanged();
	}// END: collectSettings

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	private class ListenSimulate implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frame.doExport();

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

	private class ListenGenerateXML implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			System.out.println("TODO: ListenGenerateXML");

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

}// END: class
