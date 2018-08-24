/*
 * AboutDialog.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.bss;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
@SuppressWarnings("serial")
public class AboutDialog extends JDialog {

//	+ "A BibTeX entry for LaTeX users: \n\n"
//	+ "@Article{, \n"
//	+ "\t AUTHOR = {Bielejec, Filip and Lemey, Philippe and Carvalho, Luiz and Baele, Guy and Rambaut, Andrew and Suchard, Marc A.}, \n"
//	+ "\t TITLE = {piBUSS: a parallel BEAST/BEAGLE utility for sequence simulation under complex evolutionary scenarios}, \n"
//	+ "\t JOURNAL = {BMC Bioinformatics}, \n"
//	+ "\t VOLUME = {15}, \n"
//	+ "\t YEAR = {2014}, \n" 
//	+ "\t NUMBER = {1}, \n" 
//	+ "\t PAGES = {133} \n" +
//	"}";
	
	private static final int WIDTH = 700;
	private static final int HEIGHT = 700;
	private static final int FONT_SIZE = 15;

	private static final String FILIP_BIELEJEC = "Filip Bielejec";
	private static final String ANDREW_RAMBAUT = "Andrew Rambaut";
	private static final String MARC_SUCHARD = "Marc A. Suchard";
	private static final String PHILIPPE_LEMEY = "Philippe Lemey";
	private static final String LUIZ_MAX_CARVAHLO = "Luiz Max Carvahlo";
	private static final String GUY_BAELE = "Guy Baele";

	private static final String CITATION1 = "To cite " + BeagleSequenceSimulatorApp.SHORT_NAME + " in publications, please use:";
	private static final String CITATION2 = "Bielejec, F., P. Lemey, L. Carvalho, G. Baele, A. Rambaut, and M. A. Suchard. 2014.";
	private static final String CITATION3 = "pibuss: a parallel beast/beagle utility for sequence simulation under complex evolutionary scenarios."; 
	private static final String CITATION4 = "BMC Bioinformatics 15:133";
	
	public AboutDialog() {
		initUI();
	}// END: Constructor

	public final void initUI() {

		JLabel label;
		JLabel contact;
		JLabel website;
		String addres;

		setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
		getContentPane().setBackground(Color.WHITE);
		setLocationRelativeTo(Utils.getActiveFrame());
		
		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup image
		label = new JLabel(Utils.createImageIcon(Utils.BSS_ICON));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup name
		label = new JLabel(BeagleSequenceSimulatorApp.SHORT_NAME);
		label.setFont(new Font("Serif", Font.BOLD, FONT_SIZE));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup long name
		label = new JLabel(BeagleSequenceSimulatorApp.LONG_NAME);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);

		// Setup version
		label = new JLabel("Version v" + BeagleSequenceSimulatorApp.VERSION
				+ 
//				" Prerelease" + 
				", " + BeagleSequenceSimulatorApp.DATE);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup authors
		label = new JLabel("by " + FILIP_BIELEJEC + ", " + ANDREW_RAMBAUT
				+ ", " + MARC_SUCHARD + ", " + GUY_BAELE + ", "  + LUIZ_MAX_CARVAHLO + " and "
				+ PHILIPPE_LEMEY);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));
		
		//Setup citation
		label = new JLabel(CITATION1);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);		
		
		label = new JLabel(CITATION2);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);	
		
		label = new JLabel(CITATION3);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);	
		
		label = new JLabel(CITATION4);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 2));
		label.setAlignmentX(0.5f);
		add(label);	
		
		add(Box.createRigidArea(new Dimension(0, 10)));

		// Setup about
		label = new JLabel("BEAST auxiliary software package");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		website = new JLabel();
		addres = "http://beast.community";
		website.setText("<html><p><a href=\"" + addres + "\">" + addres
				+ "</a></p></html>");
		website.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		website.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// website.setVerticalAlignment(SwingConstants.CENTER);
		// website.setHorizontalAlignment(SwingConstants.CENTER);
		website.addMouseListener(new ListenBrowse(addres));
		add(website);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Designed and developed by");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		label = new JLabel(FILIP_BIELEJEC + ", " + MARC_SUCHARD + " and " + ANDREW_RAMBAUT);
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Computational and Evolutionary Virology, KU Leuven");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		contact = new JLabel();
		addres = "filip.bielejec(AT)rega.kuleuven.be";
		contact.setText("<html><center><p><a href=\"mailto:" + addres + "\">"
				+ addres + "</a></p></center></html>");
		contact.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		contact.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// contact.setAlignmentX(0.0f);
		contact.addMouseListener(new ListenSendMail(addres));
		add(contact);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Departments of Biomathematics and Human Genetics, University of California, Los Angeles");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		contact = new JLabel();
		addres = "msuchard(AT)ucla.edu";
		contact.setText("<html><center><p><a href=\"mailto:" + addres + "\">"
				+ addres + "</a></p></center></html>");
		contact.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		contact.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// contact.setAlignmentX(0.0f);
		contact.addMouseListener(new ListenSendMail(addres));
		add(contact);

		add(Box.createRigidArea(new Dimension(0, 10)));
		
		label = new JLabel(
				"Institute of Evolutionary Biology, University of Edinburgh");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		contact = new JLabel();
		addres = "a.rambaut(AT)ed.ac.uk";
		contact.setText("<html><p><a href=\"mailto:" + addres + "\">" + addres
				+ "</a></p></html>");
		contact.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		contact.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// contact.setAlignmentX(0.5f);
		contact.addMouseListener(new ListenSendMail(addres));
		add(contact);

		add(Box.createRigidArea(new Dimension(0, 10)));

		label = new JLabel("Source code distributed under the GNU LGPL");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);

		website = new JLabel();
		addres = "http://code.google.com/p/beast-mcmc";
		website.setText("<html><p><a href=\"" + addres + "\">" + addres
				+ "</a></p></html>");
		website.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		website.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		// website.setAlignmentX(0.5f);
		website.addMouseListener(new ListenBrowse(addres));
		add(website);

		add(Box.createRigidArea(new Dimension(0, 20)));

		label = new JLabel("In case of any problems please contact your local witch doctor or a shaman.");
		label.setFont(new Font("Serif", Font.PLAIN, FONT_SIZE - 3));
		label.setAlignmentX(0.5f);
		add(label);
		
		add(Box.createRigidArea(new Dimension(0, 20)));

		JButton close = new JButton("Close");
		close.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent event) {
				dispose();
			}
		});

		close.setAlignmentX(0.5f);
		add(close);

		setModalityType(ModalityType.APPLICATION_MODAL);

		setTitle("About " + BeagleSequenceSimulatorApp.SHORT_NAME);
		setDefaultCloseOperation(DISPOSE_ON_CLOSE);
		setLocationRelativeTo(null);
		setSize(WIDTH, HEIGHT);
		setResizable(true);
	}// END: initUI

	private class ListenSendMail extends MouseAdapter {

		private String addres;

		public ListenSendMail(String addres) {
			this.addres = addres;
		}

		@Override
		public void mouseClicked(MouseEvent ev) {
			try {

				Desktop.getDesktop().mail(new URI("mailto:" + addres));

			} catch (IOException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this address in your system's standard email client.");
			} catch (URISyntaxException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this address in your system's standard email client.");
			}// END: try-catch block

		}// END: mouseClicked

	}// END: ListenSendMail

	private class ListenBrowse extends MouseAdapter {

		private String website;

		public ListenBrowse(String website) {
			this.website = website;
		}

		@Override
		public void mouseClicked(MouseEvent ev) {
			
			try {

				Desktop.getDesktop().browse(new URI(website));

			} catch (IOException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this link in your system's standard browser.");
			} catch (URISyntaxException e) {
				Utils.handleException(
						e,
						"Problem occurred while trying to open this link in your system's standard browser.");
			}// END: try-catch block

		}// END: mouseClicked

	}// END: ListenSendMail

}// END: class