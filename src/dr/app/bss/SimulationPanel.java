package dr.app.bss;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.SwingConstants;

import org.virion.jam.components.WholeNumberField;

@SuppressWarnings("serial")
public class SimulationPanel extends JPanel implements Exportable {

	 private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JLabel replicatesLabel = new JLabel("Number of sites:");
	private WholeNumberField replicatesField = new WholeNumberField(1,
			Integer.MAX_VALUE);
	
	private ButtonGroup buttonGroup;
	private JPanel panel;
	private String homogenousSimulation;
	private String heterogenousSimulation;
	private JRadioButton homogenousSimulationRadioButton;
	private JRadioButton heterogenousSimulationRadioButton;
	
	
	public SimulationPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {

		this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		// number of sites
		replicatesField.setColumns(8);
		replicatesField.setValue(data.replicateCount);
		
		optionPanel.addComponents(replicatesLabel, replicatesField);

		// simulation type
		panel = new JPanel();
		panel.setLayout(new GridLayout(2, 1));
		buttonGroup = new ButtonGroup();

		homogenousSimulation = "Homogenous simulation";
		homogenousSimulationRadioButton = new JRadioButton(homogenousSimulation);
		homogenousSimulationRadioButton.addActionListener(new ChooseSimulationTypeListener());
		homogenousSimulationRadioButton.setSelected(true);
//		frame.homogenousSimulationTypeSelected();
		buttonGroup.add(homogenousSimulationRadioButton);
		panel.add(homogenousSimulationRadioButton);
			
		heterogenousSimulation = "Heterogenous simulation";
		heterogenousSimulationRadioButton = new JRadioButton(heterogenousSimulation);
		heterogenousSimulationRadioButton.addActionListener(new ChooseSimulationTypeListener());
		heterogenousSimulationRadioButton.setSelected(true);
		buttonGroup.add(heterogenousSimulationRadioButton);
		panel.add(heterogenousSimulationRadioButton); 
		
		optionPanel.addComponentWithLabel("Choose simulation type:", panel);
		
	}// END: SimulationPanel

	public final void collectSettings() {
		data.replicateCount = replicatesField.getValue();
//		frame.fireModelChanged();
	}

	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

	public void setHomogenousSimulation() {
		homogenousSimulationRadioButton.setSelected(true);
	}// END: setHomogenousSimulation

	public void setHeterogenousSimulation() {
		heterogenousSimulationRadioButton.setSelected(true);
	}// END: setHeterogenousSimulation
	
	class ChooseSimulationTypeListener implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			if (ev.getActionCommand() == homogenousSimulation) {

				frame.homogenousSimulationTypeSelected();

			} else if (ev.getActionCommand() == heterogenousSimulation) {

				frame.heterogenousSimulationTypeSelected();

			} else {
				
				System.err.println("Unimplemented simulation type selected");
			
			}

		}// END: actionPerformed
	}// END: ChooseAnalysisTypeListener
	
	
	
}// END: class
