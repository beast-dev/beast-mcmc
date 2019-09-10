/*
 * TwoParalogGeneConversionSubstitutionModelParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evolution.datatype.PairedDataType;
import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.geneconversion.TwoParalogGeneConversionSubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class TwoParalogGeneConversionSubstitutionModelParser extends AbstractXMLObjectParser {

    private final String NAME = "twoParalogGeneConversionSubstitutionModel";
    private final String PARALOG_COUNTS = "paralogCounts";
    private final String GENE_CONVERSION_RATE = "geneConversionRate";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter paralogFrequencies = (Parameter) xo.getChild(PARALOG_COUNTS).getChild(Parameter.class);
        Parameter geneConversionRate = (Parameter) xo.getChild(GENE_CONVERSION_RATE).getChild(Parameter.class);
        BaseSubstitutionModel baseSubstitutionModel = (BaseSubstitutionModel) xo.getChild(BaseSubstitutionModel.class);
        PairedDataType dataType = (PairedDataType) xo.getChild(PairedDataType.class);
        return new TwoParalogGeneConversionSubstitutionModel(NAME + "(" + baseSubstitutionModel.getModelName() + ")",
                baseSubstitutionModel, paralogFrequencies, geneConversionRate, dataType);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PARALOG_COUNTS, Parameter.class),
            new ElementRule(GENE_CONVERSION_RATE, Parameter.class),
            new ElementRule(BaseSubstitutionModel.class),
            new ElementRule(PairedDataType.class)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TwoParalogGeneConversionSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
