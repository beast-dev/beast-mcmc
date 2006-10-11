/*
 * ConsoleFrame.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.console;

import org.virion.jam.framework.Application;
import org.virion.jam.framework.DocumentFrame;

import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.*;

public class ConsoleFrame extends DocumentFrame {
	PipedInputStream piOut;
	PipedInputStream piErr;
	PipedOutputStream poOut;
	PipedOutputStream poErr;
	JTextArea textArea = new JTextArea();

	public ConsoleFrame() throws IOException {

		super();
		// Set up System.out
		piOut = new PipedInputStream();
		poOut = new PipedOutputStream(piOut);
		System.setOut(new PrintStream(poOut, true));

		// Set up System.err
		piErr = new PipedInputStream();
		poErr = new PipedOutputStream(piErr);
		System.setErr(new PrintStream(poErr, true));

		textArea.addKeyListener(new KeyAdapter() {
			public void keyPressed(KeyEvent e) {
				if (e.getKeyCode() == KeyEvent.VK_ESCAPE ||
					(e.getKeyCode() == KeyEvent.VK_C && e.isControlDown()) ||
					(e.getKeyCode() == KeyEvent.VK_PERIOD && e.isMetaDown())) {
					escapePressed();
				}
			}
		});
	}

	protected void initializeComponents() {
		// Add a scrolling text area
		textArea.setEditable(false);
		textArea.setRows( 25 );
		textArea.setColumns( 80 );
		textArea.setEditable(false);
		textArea.setFont(new java.awt.Font("Monospaced", 0, 12));

		JScrollPane scrollPane = new JScrollPane(textArea);
		scrollPane.setVerticalScrollBarPolicy(javax.swing.JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		scrollPane.setHorizontalScrollBarPolicy(javax.swing.JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

		getContentPane().add(scrollPane, BorderLayout.CENTER);
		pack();
		setVisible(true);

		// Create reader threads
		new ReaderThread(piOut).start();
		new ReaderThread(piErr).start();
	}

	protected boolean readFromFile(File file) throws FileNotFoundException, IOException {
		throw new RuntimeException("Cannot read file");
	}

	protected boolean writeToFile(File file) throws IOException {

		textArea.write(new FileWriter(file));
		return true;
	}

    public JComponent getExportableComponent() {
		return textArea;
	}

	protected void escapePressed() {
		int option = JOptionPane.showConfirmDialog(this, "Are you sure you wish to stop?",
													"Stop",
													JOptionPane.OK_CANCEL_OPTION,
													JOptionPane.WARNING_MESSAGE);

		if (option == JOptionPane.OK_OPTION) {
			((ConsoleApplication)Application.getApplication()).doStop();
		}
	}

	public void doCopy() {
		textArea.copy();
	}

	public void doSelectAll() {
		textArea.selectAll();
	}

	class ReaderThread extends Thread {
		PipedInputStream pi;

		ReaderThread(PipedInputStream pi) {
			this.pi = pi;
		}

		public void run() {
			try {
				while (true) {
					try {
						Thread.sleep(100);
					} catch(InterruptedException ie) {}
					String input = "";
					do {
						int available = pi.available();
						if (available == 0) break;
						byte b[]=new byte[available];
						pi.read(b);
						input = input + new String(b);

					} while( !input.endsWith("\n") &&  !input.endsWith("\r\n") );

					if (input.length() > 0) {
						SwingUtilities.invokeLater(new WriteText(input));
					}
				}
			} catch (IOException e) {
			}
		}

		class WriteText implements Runnable {
			String text;

			WriteText(String text) {
				this.text = text;
			}

			public void run() {
				textArea.append(text);

				// Make sure the last line is always visible
				textArea.setCaretPosition(textArea.getDocument().getLength());

				// Keep the text area down to a certain character size
				int idealSize = 128000;
				int maxExcess = 16000 ;
				int excess = textArea.getDocument().getLength() - idealSize;
				if (excess >= maxExcess) {
					textArea.replaceRange("", 0, excess);
				}

				setDirty();
			}
		};
	}
}