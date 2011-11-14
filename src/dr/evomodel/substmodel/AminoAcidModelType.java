/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.substmodel;

import dr.util.Citation;
import jebl.evolution.substmodel.AminoAcidModel;

/**
 * @author Alexei Drummond
 */
public enum AminoAcidModelType {

    BLOSUM_62("Blosum62", "blosum62", Blosum62.INSTANCE),
    DAYHOFF("Dayhoff", "dayhoff", Dayhoff.INSTANCE),
    JTT("JTT", dr.evomodel.substmodel.JTT.INSTANCE),
    MT_REV_24("mtREV", MTREV.INSTANCE),
    CP_REV_45("cpREV", CPREV.INSTANCE),
    WAG("WAG", dr.evomodel.substmodel.WAG.INSTANCE),
    FLU("FLU", dr.evomodel.substmodel.FLU.INSTANCE);

    AminoAcidModelType(String displayName, EmpiricalRateMatrix matrix) {
        this(displayName, displayName, matrix);
    }

    AminoAcidModelType(String displayName, String xmlName, EmpiricalRateMatrix matrix) {
        this.displayName = displayName;
        this.xmlName = xmlName;
        this.matrix = matrix;
    }

    public String toString() {
        return displayName;
    }


    public String getXMLName() {
        return xmlName;
    }

    public EmpiricalRateMatrix getRateMatrixInstance() {
        return matrix;
    }

    public Citation getCitation() {
        return matrix.getCitations().get(0);
    }

    public static String[] xmlNames() {

        AminoAcidModelType[] values = values();

        String[] xmlNames = new String[values.length];
        for (int i = 0; i < values.length; i++) {
            xmlNames[i] = values[i].getXMLName();
        }
        return xmlNames;
    }

    private final String displayName, xmlName;
    private final EmpiricalRateMatrix matrix;
}
