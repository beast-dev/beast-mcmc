/*
 * ChemicalAminoAcidDecorator.java
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

import jebl.evolution.sequences.AminoAcids;

import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id: ChemicalAminoAcidDecorator.java,v 1.1 2005/11/11 16:40:41 rambaut Exp $
 */
public class ChemicalAminoAcidDecorator implements StateDecorator {
    Paint[] paints = new Paint[AminoAcids.getStateCount()];

    public ChemicalAminoAcidDecorator() {
        // Aliphatic
                paints[AminoAcids.A_STATE.getIndex()] =
                        paints[AminoAcids.V_STATE.getIndex()] =
                                paints[AminoAcids.I_STATE.getIndex()] =
                                        paints[AminoAcids.L_STATE.getIndex()] =
                                new Color(27, 4, 172);

        // Phenylalanine
        paints[AminoAcids.F_STATE.getIndex()] =
                new Color(0, 204, 255);

        // Sulphur
        paints[AminoAcids.C_STATE.getIndex()] =
                paints[AminoAcids.M_STATE.getIndex()] =
                        new Color(204, 236, 255);

        //Glycine
        paints[AminoAcids.G_STATE.getIndex()] =
                new Color(0, 255, 0);

        // Hydroxyl
        paints[AminoAcids.S_STATE.getIndex()] =
                paints[AminoAcids.T_STATE.getIndex()] =
                        new Color(137, 248, 139);

        // Tryptophan
        paints[AminoAcids.W_STATE.getIndex()] =
                new Color(204, 153, 255);

        // Tyrosine
        paints[AminoAcids.Y_STATE.getIndex()] =
                new Color(204, 255, 204);

        // Proline
        paints[AminoAcids.P_STATE.getIndex()] =
                new Color(255, 255, 0);

        // Acidic
        paints[AminoAcids.D_STATE.getIndex()] =
                paints[AminoAcids.E_STATE.getIndex()] =
                        new Color(255, 204, 0);

        // Amide
        paints[AminoAcids.N_STATE.getIndex()] =
                paints[AminoAcids.Q_STATE.getIndex()] =
                        new Color(244, 165, 4);

        // Basic
        paints[AminoAcids.H_STATE.getIndex()] =
                paints[AminoAcids.K_STATE.getIndex()] =
                        paints[AminoAcids.R_STATE.getIndex()] =
                                new Color(236, 21, 4);

        paints[AminoAcids.B_STATE.getIndex()] = Color.DARK_GRAY;
        paints[AminoAcids.Z_STATE.getIndex()] = Color.DARK_GRAY;
        paints[AminoAcids.X_STATE.getIndex()] = Color.GRAY;
        paints[AminoAcids.UNKNOWN_STATE.getIndex()] = Color.GRAY;
        paints[AminoAcids.STOP_STATE.getIndex()] = Color.GRAY;
        paints[AminoAcids.GAP_STATE.getIndex()] = Color.GRAY;
    };

    public Paint getStatePaint(int stateIndex) {
        return paints[stateIndex];
    }
}
