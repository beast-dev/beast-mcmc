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
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.text.BadLocationException;

import org.virion.jam.components.RealNumberField;

@SuppressWarnings("serial")
public class SiteRateModelPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JComboBox siteCombo;
	private RealNumberField[] siteParameterFields = new RealNumberField[BeagleSequenceSimulatorData.siteParameterNames.length];
    private JSpinner gammaCategoriesSpinner;
	
	public SiteRateModelPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) throws NumberFormatException, BadLocationException {

		this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		siteCombo = new JComboBox();
		siteCombo.setOpaque(false);

		for (String siteModel : BeagleSequenceSimulatorData.siteModels) {
			siteCombo.addItem(siteModel);
		}// END: fill loop

		siteCombo.addItemListener(new ListenSiteCombo());

		for (int i = 0; i < BeagleSequenceSimulatorData.siteParameterNames.length; i++) {
			siteParameterFields[i] = new RealNumberField();
			siteParameterFields[i].setColumns(8);
			siteParameterFields[i].setValue(data.siteParameterValues[i]);
		}// END: fill loop

		setSiteArguments();

	}// END: Constructor

	private void setSiteArguments() throws NumberFormatException, BadLocationException {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Site Rate model:"), siteCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = siteCombo.getSelectedIndex();
		
		for (int i = 0; i < data.siteParameterIndices[index].length; i++) {

			if(index == 1 && i == 0) {
				
				int k = data.siteParameterIndices[index][i];
				
				Integer initValue = Integer.valueOf(siteParameterFields[k].getText(0, 1)); 
				Integer	min = 0;
				Integer max = Integer.MAX_VALUE;
				Integer step = 1;
				
				SpinnerModel model = new SpinnerNumberModel(initValue, min, max, step);
				gammaCategoriesSpinner = new JSpinner(model);
				
				JPanel panel = new JPanel(new BorderLayout(6, 6));
				panel.add(gammaCategoriesSpinner, BorderLayout.WEST);
				panel.setOpaque(false);
				optionPanel.addComponentWithLabel(
						BeagleSequenceSimulatorData.siteParameterNames[k] + ":",
						panel);
				
			} else {
			
			int k = data.siteParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(siteParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					BeagleSequenceSimulatorData.siteParameterNames[k] + ":",
					panel);

			}// END: gama categories field check
			
		}// END: indices loop
	}// END: setSiteArguments

	private class ListenSiteCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			try {

				setSiteArguments();
				frame.fireModelChanged();

			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (BadLocationException e) {
				e.printStackTrace();
			}

		}// END: actionPerformed
	}// END: ListenSiteCombo

	public void collectSettings() {

		int index = siteCombo.getSelectedIndex();
		data.siteModel = index;
		
		for (int i = 0; i < BeagleSequenceSimulatorData.siteParameterNames.length; i++) {

			if(index == 1 && i == 0) { 
				
				data.siteParameterValues[i] = Double.valueOf(gammaCategoriesSpinner.getValue().toString()); 
						
			} else {
			
			data.siteParameterValues[i] = siteParameterFields[i].getValue();
			
			}// END: gama categories field check

		}// END: fill loop
	}// END: collectSettings
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
