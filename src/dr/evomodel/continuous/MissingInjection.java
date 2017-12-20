/*
 * MissingInjection.java
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

package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Marc A. Suchard
 */
public class MissingInjection {

    private static final String MISSING_INJECTION = "injectMissingTraits";
    private static final String MISSING_PROBABILITY = "missingProbability";
    private static final String TRAIT_NAME = dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TRAIT_NAME;

    private final List<TaxonInformation> taxonInformation;
    private final String traitName;

    public MissingInjection(Tree tree,  String traitName, double missingProbability)
            throws TaxonList.MissingAttributeException {
        taxonInformation = injectMissingValues(traitName, tree, missingProbability);
        this.traitName = traitName;
    }

    List<TaxonInformation> getTaxonInformation() {
        return taxonInformation;
    }

    public String getTraitName() {
        return traitName;
    }

    class InjectedMissingValue {
        final int index;
        final double originalValue;

        InjectedMissingValue(int index, double originalValue) {
            this.index = index;
            this.originalValue = originalValue;
        }
    }

    class TaxonInformation {
        final List<InjectedMissingValue> injectedMissingValues;
        final String newAttribute;
        int index;
        Taxon taxon;

        TaxonInformation(List<InjectedMissingValue> injectedMissingValues,
                         String newAttribute) {
            this.injectedMissingValues = injectedMissingValues;
            this.newAttribute = newAttribute;
        }
    }

    private List<TaxonInformation> injectMissingValues(String traitName,
                                                       Tree treeModel,
                                                       double missingProbability)
            throws TaxonList.MissingAttributeException {

        List<TaxonInformation> taxonInformationList = new ArrayList<TaxonInformation>();

        for (int index = 0, taxonCount = treeModel.getTaxonCount(); index < taxonCount; index++) {

            String object = (String) treeModel.getTaxonAttribute(index, traitName);
            if (object == null) {
                String taxonName = treeModel.getTaxonId(index);
                throw new TaxonList.MissingAttributeException("Trait \"" + traitName +
                        "\" not found for taxa \"" + taxonName + "\"");
            }

            TaxonInformation taxonInformation = injectForOneTaxon(object, missingProbability);
            Taxon taxon = treeModel.getTaxon(index);
            taxon.setAttribute(traitName, taxonInformation.newAttribute);

            taxonInformation.taxon = taxon;
            taxonInformation.index = index;

            taxonInformationList.add(taxonInformation);
        }

        return taxonInformationList;
    }

    private TaxonInformation injectForOneTaxon(String object, double missingProbability) {

        List<InjectedMissingValue> injectedMissingValues = new ArrayList<InjectedMissingValue>();
        StringBuilder sb = new StringBuilder();
        StringTokenizer st = new StringTokenizer(object);
        int count = st.countTokens();

        for (int j = 0; j < count; j++) {

            String oneValue = st.nextToken();
            if (!TreeTraitParserUtilities.isMissing(oneValue)) {
                if (MathUtils.nextDouble() < missingProbability) {
                    injectedMissingValues.add(new InjectedMissingValue(
                            j, new Double(oneValue)
                    ));
                    oneValue = "?";
                }
            }

            sb.append(oneValue).append(" ");
        }

        return new TaxonInformation(injectedMissingValues, sb.toString());
    }
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Tree tree = (Tree) xo.getChild(Tree.class);
            String traitName = xo.getStringAttribute(TRAIT_NAME);
            
            double missingProbability = xo.getDoubleAttribute(MISSING_PROBABILITY);
            if (missingProbability < 0.0 || missingProbability > 1.0) {
                throw new XMLParseException("Must provide a missing probability 0 <= x <= 1");
            }

            MissingInjection injector;
            try {
                injector = new MissingInjection(tree, traitName, missingProbability);
            } catch (TaxonList.MissingAttributeException e) {
                throw new XMLParseException(e.getMessage());
            }

            return injector;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return MissingInjection.class;
        }

        @Override
        public String getParserName() {
            return MISSING_INJECTION;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newStringRule(TRAIT_NAME),
                AttributeRule.newDoubleRule(MISSING_PROBABILITY),
        };
    };
}
