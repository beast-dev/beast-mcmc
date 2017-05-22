/*
 * NtdBMAParser.java
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

package dr.oldevomodelxml.substmodel;

import dr.xml.*;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.NtdBMA;
import dr.inference.model.Variable;

/**
 * @author Chieh-Hsi Wu
 * Parser for nucleotide subsitution
 */
public class NtdBMAParser extends AbstractXMLObjectParser {
    public static final String NTD_BMA = "ntdBMA";
    public static final String KAPPA = "kappa";
    public static final String TN = "tn";
    public static final String AT = "at";
    public static final String AC = "ac";
    public static final String GC = "gc";
    public static final String GT = "gt";
    public static final String MODEL_CHOOSE = "modelChoose";

    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return NTD_BMA;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {



        Variable kappa = (Variable) xo.getElementFirstChild(KAPPA);
        Variable tn = (Variable) xo.getElementFirstChild(TN);
        Variable ac = (Variable) xo.getElementFirstChild(AC);
        Variable at = (Variable) xo.getElementFirstChild(AT);
        Variable gc = (Variable) xo.getElementFirstChild(GC);
        Variable gt  = (Variable) xo.getElementFirstChild(GT);
        Variable modelChoose  = (Variable) xo.getElementFirstChild(MODEL_CHOOSE);
        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);


        return new NtdBMA(kappa, tn, ac, at, gc, gt, modelChoose, freqModel);
    }

    public Class getReturnType() {
        return NtdBMA.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(AC,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(AT,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(GC,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(GT,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(MODEL_CHOOSE,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
    };
        public String getParserDescription() {
        return "A model that allows model averaging over nucleotide substitution models.";
    }
}
