/*
 * PosteriorPredictiveXMLWriter.java
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

package dr.evomodel.treedatalikelihood.pps;

import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evoxml.AlignmentParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.SequenceParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.io.StringWriter;

/**
 * Builds the XML fragments that make up a posterior-predictive-simulation log file:
 * a wrapping root element, one-time taxa/dataType preambles, and one alignment block
 * per logged state. Has no dependency on Logger/MCMC/BEAGLE machinery, so each write*
 * method can be exercised directly against real evolution objects.
 */
public class PosteriorPredictiveXMLWriter {

    public static final String PPS = "pps";
    public static final String STATE = "state";
    public static final String PARENT_ID = "parentId";

    public String writeRootOpen() {
        StringWriter sw = new StringWriter();
        XMLWriter writer = new XMLWriter(sw);
        writer.writeOpenTag(PPS);
        writer.close();
        return sw.toString();
    }

    public String writeRootClose() {
        StringWriter sw = new StringWriter();
        XMLWriter writer = new XMLWriter(sw);
        writer.writeCloseTag(PPS);
        writer.close();
        return sw.toString();
    }

    public String writeTaxaBlock(TaxonList taxa, String id) {
        StringWriter sw = new StringWriter();
        XMLWriter writer = new XMLWriter(sw);

        writer.writeOpenTag(TaxaParser.TAXA, new Attribute.Default<String>(XMLParser.ID, id));
        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);
            writer.writeTag(TaxonParser.TAXON,
                    new Attribute.Default<String>(XMLParser.ID, taxon.getId()), true);
        }
        writer.writeCloseTag(TaxaParser.TAXA);

        writer.close();
        return sw.toString();
    }

    public String writeGeneralDataTypeBlock(GeneralDataType dataType, String id) {
        StringWriter sw = new StringWriter();
        XMLWriter writer = new XMLWriter(sw);

        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE,
                new Attribute.Default<String>(XMLParser.ID, id));
        for (int i = 0; i < dataType.getStateCount(); i++) {
            writer.writeTag(GeneralDataTypeParser.STATE,
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, dataType.getCode(i)), true);
        }
        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);

        writer.close();
        return sw.toString();
    }

    public String writeAlignmentBlock(SimpleAlignment alignment, String id, long state, String parentId) {
        StringWriter sw = new StringWriter();
        XMLWriter writer = new XMLWriter(sw);

        DataType dataType = alignment.getDataType();
        writer.writeOpenTag(AlignmentParser.ALIGNMENT, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, id),
                new Attribute.Default<String>(DataType.DATA_TYPE, dataType.getDescription()),
                new Attribute.Default<Long>(STATE, state),
                new Attribute.Default<String>(PARENT_ID, parentId)
        });

        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            writer.writeOpenTag(SequenceParser.SEQUENCE);
            writer.writeIDref(TaxonParser.TAXON, alignment.getTaxon(i).getId());
            writer.writeText(alignment.getSequence(i).getSequenceString());
            writer.writeCloseTag(SequenceParser.SEQUENCE);
        }

        writer.writeCloseTag(AlignmentParser.ALIGNMENT);

        writer.close();
        return sw.toString();
    }
}
