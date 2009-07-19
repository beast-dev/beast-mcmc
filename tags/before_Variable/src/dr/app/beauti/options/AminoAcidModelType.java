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

package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum AminoAcidModelType {

    BLOSUM_62("Blosum62", "blosum62"),
    DAYHOFF("Dayhoff", "dayhoff"),
    JTT("JTT"),
    MT_REV_24("mtREV"),
    CP_REV_45("cpREV"),
    WAG("WAG");

    AminoAcidModelType(String displayName) {
        this(displayName, displayName);
    }

    AminoAcidModelType(String displayName, String xmlName) {
        this.displayName = displayName;
        this.xmlName = xmlName;
    }

    public String toString() {
        return displayName;
    }


    public String getXMLName() {
        return xmlName;
    }

    public static String[] xmlNames() {

        AminoAcidModelType[] values = values();

        String[] xmlNames = new String[values.length];
        for (int i = 0; i < values.length; i++) {
            xmlNames[i] = values[i].getXMLName();
        }
        return xmlNames;
    }

    String displayName, xmlName;
}
