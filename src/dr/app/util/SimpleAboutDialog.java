/*
 * SimpleAboutDialog.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.util;

import jam.html.SimpleLinkListener;
import jam.util.IconUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

/**
 * An about box laid out like jam's AboutBox (the icon beside the title, the credits below) but needing no
 * jam Application behind it, closed by a click or Escape. For programs that show a dialog before, or
 * without, an application window, such as BEAST's options dialog.
 *
 * @author Andrew Rambaut
 */
public class SimpleAboutDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    /**
     * @param title      the title, as plain text
     * @param about      the credits (may be HTML)
     * @param icon       the application's icon (shown at 96 pixels), or null
     * @param background the box's colour, or null for the default
     * @param foreground the text colour, or null for the default
     */
    public SimpleAboutDialog(String title, String about, Icon icon, Color background, Color foreground) {
        // owned by the active window (the one whose menu was used), so that it opens in front of it
        super(KeyboardFocusManager.getCurrentKeyboardFocusManager().getActiveWindow(), "About", ModalityType.MODELESS);
        setDefaultCloseOperation(DISPOSE_ON_CLOSE);
        setResizable(false);

        JPanel panel = new JPanel(new BorderLayout(0, 10));
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        if (background != null) {
            panel.setBackground(background);
        }

        JLabel titleLabel = new JLabel(title, SwingConstants.CENTER);
        titleLabel.setFont(titleLabel.getFont().deriveFont(32.0f).deriveFont(Font.PLAIN));
        if (icon != null) {
            titleLabel.setIcon(IconUtils.resize(icon, 96, 96));
        }
        if (foreground != null) {
            titleLabel.setForeground(foreground);
        }
        panel.add(titleLabel, BorderLayout.NORTH);

        JEditorPane aboutPane = about.startsWith("<html>") ? new JEditorPane("text/html", about) : new JEditorPane("text/plain", about);
        aboutPane.setOpaque(false);
        aboutPane.setEditable(false);
        aboutPane.setFont(aboutPane.getFont().deriveFont(11.0f));
        if (foreground != null) {
            aboutPane.setForeground(foreground);
        }
        aboutPane.addHyperlinkListener(new SimpleLinkListener());
        panel.add(aboutPane, BorderLayout.CENTER);

        MouseAdapter closer = new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                dispose();
            }
        };
        panel.addMouseListener(closer);
        aboutPane.addMouseListener(closer);
        getRootPane().registerKeyboardAction(new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                dispose();
            }
        }, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_IN_FOCUSED_WINDOW);

        setContentPane(panel);
        pack();
        setLocationRelativeTo(getOwner());
    }

    @Override
    public void setVisible(boolean visible) {
        super.setVisible(visible);
        if (visible) {
            toFront();
        }
    }
}
