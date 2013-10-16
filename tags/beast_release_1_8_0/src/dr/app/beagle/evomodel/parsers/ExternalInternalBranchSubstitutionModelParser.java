/*
 * ExternalInternalBranchSubstitutionModelParser.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.ExternalInternalBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

@Deprecated

public class ExternalInternalBranchSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "tipBranchSubstitutionModel";
    public static final String INTERNAL = "internal";
    public static final String EXTERNAL = "external";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<SubstitutionModel> modelList = new ArrayList<SubstitutionModel>();
        modelList.add((SubstitutionModel) xo.getElementFirstChild(INTERNAL));
        modelList.add((SubstitutionModel) xo.getElementFirstChild(EXTERNAL));

        List<FrequencyModel> freqList = new ArrayList<FrequencyModel>();
        freqList.add((FrequencyModel) xo.getChild(FrequencyModel.class));
        return new ExternalInternalBranchSubstitutionModel(modelList, freqList);
    }

    /**
     * @return an array of syntax rules required by this element.
     *         Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(INTERNAL, SubstitutionModel.class),
                new ElementRule(EXTERNAL, SubstitutionModel.class),               
                new ElementRule(FrequencyModel.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "A branch site model that uses different site models on the external and internal branches";
    }

    @Override
    public Class getReturnType() {
        return BranchSubstitutionModel.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    public String getParserName() {
        return PARSER_NAME;
    }
}
