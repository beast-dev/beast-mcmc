/*
 * DuplicationBranchSubstitutionModel.java
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

import dr.evomodel.substmodel.geneconversion.DuplicationBranchSubstitutionModel;
import dr.evomodel.substmodel.geneconversion.PairedParalogGeneConversionSubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class DuplicationBranchSubstitutionModelParser extends AbstractXMLObjectParser{

    private final String NAME = "duplicationBranchSubstitutionModel";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter timeUntilDuplicationProportion = (Parameter) xo.getChild(Parameter.class);
        List<PairedParalogGeneConversionSubstitutionModel> substitutionModels = new ArrayList<>();
        for (PairedParalogGeneConversionSubstitutionModel substitutionModel : xo.getAllChildren(PairedParalogGeneConversionSubstitutionModel.class)) {
            substitutionModels.add(substitutionModel);
        }
        return new DuplicationBranchSubstitutionModel(NAME, substitutionModels, timeUntilDuplicationProportion);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(PairedParalogGeneConversionSubstitutionModel.class, 2, 2),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return DuplicationBranchSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
