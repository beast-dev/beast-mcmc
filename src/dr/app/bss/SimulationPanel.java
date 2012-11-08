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
			printDataList();
			
		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates
	
	private void printDataList() {

		int row = 1;
		for (PartitionData data : dataList) {

			System.out.println("Partition: " + row);
			System.out.println("\tReplications: " + dataList.sequenceLength);
			System.out.println("\tFrom: " + data.from);
			System.out.println("\tTo: " + data.to);
			System.out.println("\tEvery: " + data.every);
			System.out.println("\tTree model: " + data.treeFile.getName());
			System.out.println("\tSubstitution model: " + PartitionData.substitutionModels[data.substitutionModel]);
			System.out.println("\tSite rate model: " + PartitionData.siteModels[data.siteModel]);
			System.out.println("\tClock rate model: " + PartitionData.clockModels[data.clockModel]);
			System.out.println("\tFrequency model: " + PartitionData.frequencyModels[data.frequencyModel]);

			System.out.println("Possible trees: ");
			for (int i = 0; i < dataList.treeFilesList.size(); i++) {
				System.out.println(dataList.treeFilesList.get(i).getName());
			}
			
			row++;
		}// END: data list loop

	}// END: printDataList
	
}// END: class
