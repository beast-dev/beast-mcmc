/*
 * MutationDeathModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.MutationDeathModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.datatype.MutationDeathType;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexander V. Alekseyenko
 * @author Marc A. Suchard
 */
public class MutationDeathModelParser extends AbstractXMLObjectParser {

    public static final String MD_MODEL = "mutationDeathModel";
    public static final String MUTATION_RATE = "mutationRate";

    public String getParserName() {
        return MD_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter dummyFreqParameter;
        Parameter delParam = (Parameter) xo.getChild(Parameter.class);

        Logger.getLogger("dr.evomodel").info("Creating MutationDeath substitution model.\n\tInitial death rate is "
                + delParam.getParameterValue(0));

        MutationDeathType dT = (MutationDeathType) xo.getChild(MutationDeathType.class);

        SubstitutionModel evoModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        if (evoModel == null) {  // Assuming pure survival model
            Logger.getLogger("dr.evomodel").info("\tSubstitutionModel not provided assuming pure death/survival model.");
            dummyFreqParameter = new Parameter.Default(new double[]{1.0, 0.0});
        } else {
            dummyFreqParameter = new Parameter.Default(dT.getStateCount());
            double freqs[] = evoModel.getFrequencyModel().getFrequencies();
            for (int i = 0; i < freqs.length; ++i) {
                dummyFreqParameter.setParameterValueQuietly(i, freqs[i]);
            }
            dummyFreqParameter.setParameterValueQuietly(dT.getStateCount() - 1, 0.0);
        }

        FrequencyModel dummyFrequencies = new FrequencyModel(dT, dummyFreqParameter);

        Parameter mutationRate;

        if (xo.hasChildNamed(MUTATION_RATE)) {
            mutationRate = (Parameter) xo.getElementFirstChild(MUTATION_RATE);
        } else {
            mutationRate = new Parameter.Default(new double[]{1.0});
        }
        Logger.getLogger("dr.evomodel").info("\tInitial mutation rate is " + mutationRate.getParameterValue(0));

        return new MutationDeathModel(delParam, dT, evoModel, dummyFrequencies, mutationRate);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the MutationDeath model of CTMC evolution with deletions.";
    }

    public Class getReturnType() {
        return MutationDeathModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(Parameter.class),
            new ElementRule(MutationDeathType.class),
            new ElementRule(MUTATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };

}
