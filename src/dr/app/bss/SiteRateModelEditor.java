package dr.app.bss;

import jam.panels.OptionsPanel;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.text.BadLocationException;

import org.virion.jam.components.RealNumberField;

public class SiteRateModelEditor {

	// Data
	private PartitionDataList dataList = null;
	private int row;
	
	// Settings	
	private OptionsPanel optionPanel;
	private JComboBox siteCombo;
	private RealNumberField[] siteParameterFields = new RealNumberField[PartitionData.siteParameterNames.length];
    private JSpinner gammaCategoriesSpinner;
	
	//Buttons
	private JButton done;
	private JButton cancel;
	
	// Window
	private JDialog window;
	private Frame owner;
    
	public SiteRateModelEditor(PartitionDataList dataList, int row) throws NumberFormatException, BadLocationException {

		this.dataList = dataList;
		this.row = row;
		
		window = new JDialog(owner, "Setup site rate model for partition " + (row + 1));
		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);

		siteCombo = new JComboBox();
		siteCombo.setOpaque(false);

		for (String siteModel : PartitionData.siteModels) {
			siteCombo.addItem(siteModel);
		}// END: fill loop

		siteCombo.addItemListener(new ListenSiteCombo());

		for (int i = 0; i < PartitionData.siteParameterNames.length; i++) {
			siteParameterFields[i] = new RealNumberField();
			siteParameterFields[i].setColumns(8);
			siteParameterFields[i].setValue(dataList.get(0).siteParameterValues[i]);
		}// END: fill loop

		setSiteArguments();

		// Buttons
		JPanel buttonsHolder = new JPanel();
		buttonsHolder.setOpaque(false);
		
		cancel = new JButton("Cancel", BeagleSequenceSimulatorApp.closeIcon);
		cancel.addActionListener(new ListenCancel());
		buttonsHolder.add(cancel);
		
		done = new JButton("Done", BeagleSequenceSimulatorApp.doneIcon);
		done.addActionListener(new ListenOk());
		buttonsHolder.add(done);
		
		// Window
		owner = Utils.getActiveFrame();
		window.setLocationRelativeTo(owner);
		window.getContentPane().setLayout(new BorderLayout());
		window.getContentPane().add(optionPanel, BorderLayout.CENTER);
		window.getContentPane().add(buttonsHolder, BorderLayout.SOUTH);
		window.pack();
	}// END: Constructor

	private void setSiteArguments() throws NumberFormatException, BadLocationException {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Site Rate model:"), siteCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = siteCombo.getSelectedIndex();
		
		for (int i = 0; i < dataList.get(0).siteParameterIndices[index].length; i++) {

			if(index == 1 && i == 0) {
				
				int k = dataList.get(0).siteParameterIndices[index][i];
				
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
						PartitionData.siteParameterNames[k] + ":",
						panel);
				
			} else {
			
			int k = dataList.get(0).siteParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(siteParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					PartitionData.siteParameterNames[k] + ":",
					panel);

			}// END: gama categories field check
			
		}// END: indices loop
		
		window.validate();
		window.repaint();
	}// END: setSiteArguments

	private class ListenSiteCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			try {

				setSiteArguments();

			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (BadLocationException e) {
				e.printStackTrace();
			}

		}// END: actionPerformed
	}// END: ListenSiteCombo

	public void collectSettings() {

		int index = siteCombo.getSelectedIndex();
		dataList.get(row).siteModel = index;
		
		for (int i = 0; i < PartitionData.siteParameterNames.length; i++) {

			if(index == 1 && i == 0) { 
				
				dataList.get(row).siteParameterValues[i] = Double.valueOf(gammaCategoriesSpinner.getValue().toString()); 
						
			} else {
			
				dataList.get(0).siteParameterValues[i] = siteParameterFields[i].getValue();
			
			}// END: gama categories field check

		}// END: fill loop
	}// END: collectSettings
	
	private class ListenOk implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			collectSettings();
			
		}// END: actionPerformed
	}// END: ListenSaveLocationCoordinates
	
	private class ListenCancel implements ActionListener {
		public void actionPerformed(ActionEvent ev) {

			window.setVisible(false);
			
		}// END: actionPerformed
	}// END: ListenCancel
	
	public void launch() {
		window.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		window.setSize(new Dimension(450, 400));
		window.setMinimumSize(new Dimension(100, 100));
		window.setResizable(true);
		window.setVisible(true);
	}//END: launch

}// END: class
