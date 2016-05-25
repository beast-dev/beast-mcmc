/*
 * GTOPO30Browser.java
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

package dr.geo;

import dr.app.gui.ColorFunction;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.io.IOException;

/**
 * @author Alexei Drummond
 */
public class GTOPO30Browser extends JFrame implements ChangeListener {

    GTOPO30Panel gtopo30panel;
    JScrollPane scrollPane;

    public GTOPO30Browser(String[] tilefiles) throws IOException {
        super("GTOPO30 Browser");

        ColorFunction function = new ColorFunction(
                new Color[]{Color.blue, Color.yellow, Color.green.darker(), Color.orange.darker(), Color.white, Color.pink},
                new float[]{-410, 0, 100, 1500, 4000, 8800});

        gtopo30panel = new GTOPO30Panel(tilefiles, function);


        scrollPane = new JScrollPane(gtopo30panel);

        getContentPane().add(BorderLayout.CENTER, scrollPane);

        JSlider scale = new JSlider(JSlider.HORIZONTAL, 1, 100, 10);
        scale.addChangeListener(this);

        getContentPane().add(BorderLayout.SOUTH, scale);
    }

    public static void main(String[] args) throws IOException {

        GTOPO30Browser browser = new GTOPO30Browser(args);
        browser.setSize(800, 800);
        browser.setVisible(true);
    }


    public void stateChanged(ChangeEvent e) {
        JSlider source = (JSlider) e.getSource();
        if (!source.getValueIsAdjusting()) {
            int scale = source.getValue();
            gtopo30panel.setScale((double) scale / 100.0f);
            scrollPane.repaint();
        }
    }
}
