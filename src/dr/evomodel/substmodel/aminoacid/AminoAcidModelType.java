/*
 * AminoAcidModelType.java
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

package dr.evomodel.substmodel.aminoacid;

import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public enum AminoAcidModelType {

    BLOSUM_62("Blosum62", "blosum62", Blosum62.INSTANCE),
    DAYHOFF("Dayhoff", "dayhoff", Dayhoff.INSTANCE),
    JTT("JTT", dr.evomodel.substmodel.aminoacid.JTT.INSTANCE),
    MT_REV_24("mtREV", MTREV.INSTANCE),
    CP_REV_45("cpREV", CPREV.INSTANCE),
    WAG("WAG", dr.evomodel.substmodel.aminoacid.WAG.INSTANCE),
    LG("LG", dr.evomodel.substmodel.aminoacid.LG.INSTANCE),
    FLU("FLU", dr.evomodel.substmodel.aminoacid.FLU.INSTANCE),
    MTVER("mtVer", dr.evomodel.substmodel.aminoacid.MTVER.INSTANCE),
    MTPRO("mtPro", dr.evomodel.substmodel.aminoacid.MTPRO.INSTANCE),
    MTMET("mtMet", dr.evomodel.substmodel.aminoacid.MTMET.INSTANCE),
    MTINV("mtInv", dr.evomodel.substmodel.aminoacid.MTINV.INSTANCE),
    MTDEU("mtDeu", dr.evomodel.substmodel.aminoacid.MTDEU.INSTANCE),
    MTMAM("mtMam", dr.evomodel.substmodel.aminoacid.MTMAM.INSTANCE),
    Q_3DI_AF("Q.3Di.AF", ThreeDiAF.INSTANCE),
    Q_3DI_LLM("Q.3Di.LLM", ThreeDiLLM.INSTANCE),
    Q_3DI_FS_LOGODDS("Q.3Di.FSlogodds", ThreeDiFSlogodds.INSTANCE, "ThreeDi");

    AminoAcidModelType(String displayName, EmpiricalRateMatrix matrix, String... legacyXmlNames) {
        this(displayName, displayName, matrix, legacyXmlNames);
    }

    AminoAcidModelType(String displayName, String xmlName, EmpiricalRateMatrix matrix, String... legacyXmlNames) {
        this.displayName = displayName;
        this.xmlName = xmlName;
        this.matrix = matrix;
        this.legacyXmlNames = legacyXmlNames;
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

    /**
     * @return the model whose current or legacy XML name matches, or null if there is none.
     */
    public static AminoAcidModelType fromXMLName(String xmlName) {

        for (AminoAcidModelType type : values()) {
            if (type.xmlName.equals(xmlName)) {
                return type;
            }
            for (String legacyXmlName : type.legacyXmlNames) {
                if (legacyXmlName.equals(xmlName)) {
                    return type;
                }
            }
        }
        return null;
    }

    /**
     * @return every XML name accepted for a model, legacy names of renamed models included.
     * This is what the XML syntax rules validate against, so a name omitted here is rejected
     * before the parser is ever reached.
     */
    public static String[] xmlNames() {

        List<String> xmlNames = new ArrayList<String>();
        for (AminoAcidModelType type : values()) {
            xmlNames.add(type.getXMLName());
            Collections.addAll(xmlNames, type.legacyXmlNames);
        }
        return xmlNames.toArray(new String[0]);
    }

    private final String displayName, xmlName;
    private final String[] legacyXmlNames;
    private final EmpiricalRateMatrix matrix;
}
