/*
 * AnalysisMenuFactory.java
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

package dr.app.tracer;

import org.virion.jam.framework.MenuFactory;
import org.virion.jam.framework.AbstractFrame;

import javax.swing.*;

/**
 * @author rambaut
 *         Date: Feb 24, 2005
 *         Time: 5:12:11 PM
 */
public class AnalysisMenuFactory implements MenuFactory {

    public static final String DEMOGRAPHIC_RECONSTRUCTION = "Demographic Reconstruction...";
    public static final String BAYESIAN_SKYLINE_RECONSTRUCTION = "Bayesian Skyline Reconstruction...";

    public String getMenuName() {
        return "Analysis";
    }

    public void populateMenu(JMenu menu, AbstractFrame frame) {
        JMenuItem item;

        if (frame instanceof AnalysisMenuHandler) {
            item = new JMenuItem(((AnalysisMenuHandler)frame).getDemographicAction());
            menu.add(item);

            item = new JMenuItem(((AnalysisMenuHandler)frame).getBayesianSkylineAction());
            menu.add(item);
        } else {
            item = new JMenuItem(DEMOGRAPHIC_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);

            item = new JMenuItem(BAYESIAN_SKYLINE_RECONSTRUCTION);
            item.setEnabled(false);
            menu.add(item);
        }

    }

    public int getPreferredAlignment() {
        return LEFT;
    }
}
