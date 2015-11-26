/*
 * GTOPO30RandomWalk.java
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
import java.awt.*;
import java.io.IOException;
import java.util.Random;

/**
 * @author Alexei Drummond
 */
public class GTOPO30RandomWalk {

    public static void main(String[] tilefiles) throws IOException {

        Color[] colors = {Color.blue, Color.yellow, Color.green.darker(), Color.orange.darker(), Color.white, Color.pink};

        ColorFunction function = new ColorFunction(
                colors,
                new float[]{-410, 0, 100, 1500, 4000, 8800});

        GTOPO30Panel gtopo30panel = new GTOPO30Panel(tilefiles, function);

        Random random = new Random();

        RateMatrix rates = gtopo30panel.getRates(100);

        gtopo30panel.setSize(gtopo30panel.latticeWidth(), gtopo30panel.latticeHeight());
        gtopo30panel.layout();

        InhomogeneousRandomWalk walk = new InhomogeneousRandomWalk(gtopo30panel, new Location(0, 0), random, rates);

        Location start = gtopo30panel.getLocation(41.9, 12.5);
        System.out.println("Start = " + start);

        for (int j = 0; j < 1000; j++) {
            walk.simulate(start, 25000);
        }
        JFrame frame = new JFrame();

        JScrollPane scrollPane = new JScrollPane(walk);

        frame.getContentPane().add(BorderLayout.CENTER, scrollPane);
        frame.setSize(1800, 1100);
        frame.setVisible(true);
    }
}
