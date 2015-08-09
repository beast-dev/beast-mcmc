/*
 * NucleotideDecorator.java
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

package dr.app.beauti.alignmentviewer;

import jebl.evolution.sequences.Nucleotides;

import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id: NucleotideDecorator.java,v 1.1 2005/11/11 16:40:41 rambaut Exp $
 */
public class NucleotideDecorator implements StateDecorator {

	public static final NucleotideDecorator INSTANCE = new NucleotideDecorator();

    Paint[] paints = new Paint[Nucleotides.getStateCount()];

    public NucleotideDecorator() {
        paints[Nucleotides.A_STATE.getIndex()] = Color.RED;
        paints[Nucleotides.C_STATE.getIndex()] = Color.BLUE;
        paints[Nucleotides.G_STATE.getIndex()] = Color.BLACK;
        paints[Nucleotides.T_STATE.getIndex()] = Color.GREEN;
        paints[Nucleotides.R_STATE.getIndex()] = Color.MAGENTA;
        paints[Nucleotides.Y_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.M_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.W_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.S_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.K_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.B_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.D_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.H_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.V_STATE.getIndex()] = Color.DARK_GRAY;
        paints[Nucleotides.N_STATE.getIndex()] = Color.GRAY;
        paints[Nucleotides.UNKNOWN_STATE.getIndex()] = Color.GRAY;
        paints[Nucleotides.GAP_STATE.getIndex()] = Color.GRAY;
    }

    public Paint getStatePaint(int stateIndex) {
        if (stateIndex >= paints.length) {
            return paints[paints.length - 1];
        }
        return paints[stateIndex];
    }
}
