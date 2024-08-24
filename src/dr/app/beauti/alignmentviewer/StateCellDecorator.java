/*
 * StateCellDecorator.java
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

package dr.app.beauti.alignmentviewer;

import java.awt.*;

/**
 * @author Andrew Rambaut
 */
public class StateCellDecorator implements CellDecorator {

    public StateCellDecorator(StateDecorator stateDecorator, boolean inverted) {
        if (inverted) {
            foregroundDecorator = new StateDecorator() {
                public Paint getStatePaint(int stateIndex) { return Color.WHITE; }
            };
            backgroundDecorator = stateDecorator;
        } else {
            foregroundDecorator = stateDecorator;
            backgroundDecorator = new StateDecorator() {
                public Paint getStatePaint(int stateIndex) { return Color.WHITE; }
            };

        }
    }

    public Paint getCellForeground(int row, int column, int state) {
        return foregroundDecorator.getStatePaint(state);
    }

    public Paint getCellBackground(int row, int column, int state) {
        return backgroundDecorator.getStatePaint(state);
    }

    private final StateDecorator foregroundDecorator;
    private final StateDecorator backgroundDecorator;

}