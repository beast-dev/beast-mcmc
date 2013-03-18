package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.SwingConstants;

import dr.app.gui.components.WholeNumberField;


@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

	public static final int FIRST_SIMULATION_TYPE = 0;
	public static final int SECOND_SIMULATION_TYPE = 1;
	public int simulationType;
	
	private MainFrame frame;
	private PartitionDataList dataList;
	private OptionsPanel optionPanel;

	private WholeNumberField siteCountField;

	private WholeNumberField simulationsNumberField;

	// Buttons
	private JButton simulate;
	private JButton generateXML;

	// Radio buttons
	private JRadioButton numberOfSimulationsRadioButton;
	private JRadioButton simulateForEachTreeRadioButton;

	// Strings for radio buttons
	private String firstSimulationType;
	private String secondSimulationType;

	public SimulationPanel(final MainFrame frame,
			final PartitionDataList dataList) {

		this.frame = frame;
		this.dataList = dataList;

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		// number of sites
		siteCountField = new WholeNumberField(1, Integer.MAX_VALUE);
		siteCountField.setColumns(8);
		siteCountField.setValue(dataList.siteCount);
		optionPanel.addComponentWithLabel("Number of sites:",
				siteCountField);

		// Simulation Type
		JPanel simulationTypeHolder = new JPanel();
		simulationTypeHolder.setOpaque(false);
		GridBagLayout layout = new GridBagLayout();
		GridBagConstraints constraints = new GridBagConstraints();
		constraints.anchor = GridBagConstraints.BASELINE_LEADING;
		simulationTypeHolder.setLayout(layout);

		// Add first radio button
		firstSimulationType = new String("Number of simulations");
		numberOfSimulationsRadioButton = new JRadioButton();
		numberOfSimulationsRadioButton.setActionCommand(firstSimulationType);
		numberOfSimulationsRadioButton
				.setToolTipText("Generate specified number of datasets.");
		numberOfSimulationsRadioButton
				.addActionListener(new ChooseAnalysisTypeListener());
		constraints.gridx = 0;
		constraints.gridy = 0;
		simulationTypeHolder.add(numberOfSimulationsRadioButton, constraints);

		// Add label
		JLabel simulationsNumberLabel = new JLabel(firstSimulationType + ":");
		constraints.gridx = 1;
		simulationTypeHolder.add(simulationsNumberLabel, constraints);

		// Add number field
		simulationsNumberField = new WholeNumberField(1, Integer.MAX_VALUE);
		simulationsNumberField.setColumns(8);
		simulationsNumberField.setValue(dataList.simulationsCount);
		constraints.gridx = 2;
		simulationTypeHolder.add(simulationsNumberField, constraints);

		// Add second radio button
		secondSimulationType = new String("Simulate for each tree");
		simulateForEachTreeRadioButton = new JRadioButton();
		simulateForEachTreeRadioButton.setActionCommand(secondSimulationType);
		simulateForEachTreeRadioButton
				.setToolTipText("Generate one dataset for each tree topology in .trees file.");
		simulateForEachTreeRadioButton
				.addActionListener(new ChooseAnalysisTypeListener());
		constraints.gridx = 0;
		constraints.gridy = 1;
		simulationTypeHolder.add(simulateForEachTreeRadioButton, constraints);

		// Add label
		JLabel simulateForEachTreeLabel = new JLabel(secondSimulationType);
		constraints.gridx = 1;
		simulationTypeHolder.add(simulateForEachTreeLabel, constraints);

		// Group radio buttons
		ButtonGroup buttonGroup = new ButtonGroup();
		buttonGroup.add(numberOfSimulationsRadioButton);
		buttonGroup.add(simulateForEachTreeRadioButton);

		optionPanel.addComponentWithLabel("Simulation type:",
				simulationTypeHolder);

		// Buttons holder
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);

		// simulate button
		simulate = new JButton("Simulate",
				BeagleSequenceSimulatorApp.biohazardIcon);
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

		setFirstSimulationType();

	}// END: SimulationPanel

	public final void collectSettings() {

		dataList.siteCount = siteCountField.getValue();
		dataList.simulationsCount = simulationsNumberField.getValue();

		// frame.fireModelChanged();
	}// END: collectSettings

	private class ListenSimulate implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frame.doExport();

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

	private class ListenGenerateXML implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			frame.doGenerateXML();

		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates

	private void setFirstSimulationType() {
		
		numberOfSimulationsRadioButton.setSelected(true);
		
		simulationsNumberField.setEnabled(true);
		generateXML.setEnabled(true);
		simulationType = FIRST_SIMULATION_TYPE;

		frame.disableTreesFileButton();
		frame.enableTreeFileButton();
		frame.showTreeColumn();
	}// END: setFirstSimulationType

	private void setSecondSimulationType() {
		simulationsNumberField.setEnabled(false);
		generateXML.setEnabled(false);
		simulationType = SECOND_SIMULATION_TYPE;

		frame.enableTreesFileButton();
		frame.disableTreeFileButton();
		frame.hideTreeColumn();
	}// END: setSecondSimulationType

	public void setBusy() {
		simulate.setEnabled(false);
		generateXML.setEnabled(false);
	}// END: setBusy

	public void setIdle() {
		simulate.setEnabled(true);
		generateXML.setEnabled(true);
	}// END: setIdle
	
	class ChooseAnalysisTypeListener implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			if (ev.getActionCommand() == firstSimulationType) {

				setFirstSimulationType();
				frame.setStatus(firstSimulationType + " selected");

			} else if (ev.getActionCommand() == secondSimulationType) {

				setSecondSimulationType();
				frame.setStatus(secondSimulationType + " selected");

			} else {
				throw new RuntimeException(
						"Unimplemented analysis type selected");
			}

		}// END: actionPerformed
	}// END: ChooseAnalysisTypeListener

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	public void setDataList(PartitionDataList  dataList) {
		this.dataList = dataList;
	}
	
}// END: class
