/*
 * XMLMerger.java
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

package dr.app.bfe;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Package: XMLMerger
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Apr 15, 2009
 *         Time: 2:42:38 AM
 */
public class XMLMerger implements ActionListener, ListSelectionListener {

    private JButton linkButton;
    private JPanel panel;
    private XMLViewer XMLViewer1;
    private XMLViewer XMLViewer2;
    private JList linkingList;
    private JButton unlinkButton;
    private JTextField linkingField;
    private JFrame frame;

    public XMLMerger() {

    }

    public XMLMerger(boolean standalone) {
        XMLViewer1.getIdElementList().addListSelectionListener(this);
        if (standalone) {
            frame = new JFrame("XML Structure Viewer");
            frame.setContentPane(panel);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.pack();
            frame.setVisible(true);
        }
    }

    public static void main(String[] args) {
        XMLMerger z = new XMLMerger(true);
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == linkButton) {
            String id1 = (String) XMLViewer1.getIdElementList().getSelectedValue();
            String id2 = (String) XMLViewer2.getIdElementList().getSelectedValue();
            String linkName = (String) linkingField.getSelectedText().trim();
            if (linkName.equals("")) {  // linking name must be non-empty
                linkName = "" + id1;
            }
            if (id1 != null && id2 != null) {
                //TODO: need to check if the id are already linked to something else
                linkingList.updateUI();
            }
        } else if (e.getSource() == unlinkButton) {

        }
    }

    public void valueChanged(ListSelectionEvent listSelectionEvent) {
        if (listSelectionEvent.getSource() == linkingList) {
            linkButton.setEnabled(false);
            unlinkButton.setEnabled(true);
        } else if (listSelectionEvent.getSource() == XMLViewer1.getIdElementList() || listSelectionEvent.getSource() == XMLViewer2.getIdElementList()) {
            if (XMLViewer1.getIdElementList().getSelectedIndex() != -1 && XMLViewer2.getIdElementList().getSelectedIndex() != -1) {
                linkButton.setEnabled(true);
            }

        }
    }
}
