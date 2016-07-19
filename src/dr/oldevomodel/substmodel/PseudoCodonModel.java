/*
 * PseudoCodonModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * PseudoCodonModel - provides an approximation to the codon substitution model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: PseudoCodonModel.java,v 1.9 2005/05/24 20:25:58 rambaut Exp $
 */

public class PseudoCodonModel extends AbstractModel {

    //
    // Public stuff
    //

    public static final String PSEUDO_CODON_MODEL = "pseudoCodonModel";
    public static final String MU = "mu";
    public static final String OMEGA = "omega";
    public static final String KAPPA = "kappa";
    public static final String FIRST_POSITION = "firstPosition";
    public static final String SECOND_POSITION = "secondPosition";
    public static final String THIRD_POSITION = "thirdPosition";

    /**
     * Constructor
     */
    public PseudoCodonModel(GammaSiteModel siteModel1,
                            GammaSiteModel siteModel2,
                            GammaSiteModel siteModel3,
                            Parameter muParameter,
                            Parameter omegaParameter,
                            Parameter kappaParameter,
                            FrequencyModel freqModel) {

        super(PSEUDO_CODON_MODEL);

        this.gtr1 = (GTR) siteModel1.getSubstitutionModel();
        this.siteModel1 = siteModel1;

        this.gtr2 = (GTR) siteModel2.getSubstitutionModel();
        this.siteModel2 = siteModel2;

        this.gtr3 = (GTR) siteModel3.getSubstitutionModel();
        this.siteModel3 = siteModel3;

        if (freqModel.getDataType() != Nucleotides.INSTANCE) {
            throw new IllegalArgumentException("Datatypes do not match!");
        }

        this.frequencyModel = freqModel;
        addModel(frequencyModel);

        this.muParameter = muParameter;
        muParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        addVariable(muParameter);

        this.omegaParameter = omegaParameter;
        omegaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        addVariable(omegaParameter);

        this.kappaParameter = kappaParameter;
        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        addVariable(kappaParameter);

        excludeStopCodons = true;

        calculateSubstitutionModel();
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // frequencyModel changed!
        calculateSubstitutionModel();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // parameter changed!
        calculateSubstitutionModel();
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    //
    // Private stuff
    //

    private void calculateSubstitutionModel() {
        double omega = omegaParameter.getParameterValue(0);
        double kappa = kappaParameter.getParameterValue(0);
        double mu = muParameter.getParameterValue(0);

        double[] frequencies = frequencyModel.getFrequencies();
        double piAC = frequencies[0] * frequencies[1];
        double piAG = frequencies[0] * frequencies[2];
        double piAT = frequencies[0] * frequencies[3];
        double piCG = frequencies[1] * frequencies[2];
        double piCT = frequencies[1] * frequencies[3];
        double piGT = frequencies[2] * frequencies[3];

        double[][] rates = new double[3][6];

        if (excludeStopCodons) {
            // A-C
            rates[0][0] = ((omega * 14.0 / 16.0) + (2.0 / 16.0));
            rates[1][0] = omega;
            rates[2][0] = ((omega * 5.0 / 14.0) + (9.0 / 14.0));

            // A-G
            rates[0][1] = omega * kappa;
            rates[1][1] = rates[0][1];
            rates[2][1] = ((omega / 14.0) + (13.0 / 14.0)) * kappa;

            // A-T
            rates[0][2] = omega;
            rates[1][2] = rates[0][2];
            rates[2][2] = ((omega * 5.0 / 14.0) + (9.0 / 14.0));

            // C-G
            rates[0][3] = omega;
            rates[1][3] = rates[0][3];
            rates[2][3] = ((omega * 7.0 / 15.0) + (8.0 / 15.0));

            // C-T
            rates[0][4] = ((omega * 11.0 / 13.0) + (2.0 / 13.0)) * kappa;
            rates[1][4] = omega * kappa;
            rates[2][4] = kappa;

            // G-T
            rates[0][5] = omega;
            rates[1][5] = rates[0][5];
            rates[2][5] = ((omega * 7.0 / 15.0) + (8.0 / 15.0));
        } else {
            // A-C
            rates[0][0] = ((omega * 14.0 / 16.0) + (2.0 / 16.0));
            rates[1][0] = omega;
            rates[2][0] = ((omega * 5.0 / 16.0) + (9.0 / 16.0));

            // A-G
            rates[0][1] = omega * kappa;
            rates[1][1] = rates[0][1];
            rates[2][1] = ((omega / 16.0) + (13.0 / 16.0)) * kappa;

            // A-T
            rates[0][2] = (omega * 13.0 / 16.0);
            rates[1][2] = omega;
            rates[2][2] = ((omega * 5.0 / 16.0) + (9.0 / 16.0));

            // C-G
            rates[0][3] = omega;
            rates[1][3] = rates[0][3];
            rates[2][3] = ((omega * 7.0 / 16.0) + (8.0 / 16.0));

            // C-T
            rates[0][4] = ((omega * 11.0 / 16.0) + (2.0 / 16.0)) * kappa;
            rates[1][4] = omega * kappa;
            rates[2][4] = kappa;

            // G-T
            rates[0][5] = (omega * 13.0 / 16.0);
            rates[1][5] = omega;
            rates[2][5] = ((omega * 7.0 / 16.0) + (8.0 / 16.0));
        }

        double sumRates1 = (rates[0][0] * piAC) +
                (rates[0][1] * piAG) +
                (rates[0][2] * piAT) +
                (rates[0][3] * piCG) +
                (rates[0][4] * piCT) +
                (rates[0][5] * piGT);
        double sumRates2 = (rates[1][0] * piAC) +
                (rates[1][1] * piAG) +
                (rates[1][2] * piAT) +
                (rates[1][3] * piCG) +
                (rates[1][4] * piCT) +
                (rates[1][5] * piGT);
        double sumRates3 = (rates[2][0] * piAC) +
                (rates[2][1] * piAG) +
                (rates[2][2] * piAT) +
                (rates[2][3] * piCG) +
                (rates[2][4] * piCT) +
                (rates[2][5] * piGT);

        double f = (mu * 3.0) / (sumRates1 + sumRates2 + sumRates3);

        double mu1 = f * sumRates1;
        double mu2 = f * sumRates2;
        double mu3 = f * sumRates3;

        gtr1.setAbsoluteRates(rates[0], 4);
        siteModel1.setMu(mu1);
        gtr2.setAbsoluteRates(rates[1], 4);
        siteModel2.setMu(mu2);
        gtr3.setAbsoluteRates(rates[2], 4);
        siteModel3.setMu(mu3);
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    public static XMLObjectParser PSEUDO_CODON_MODEL_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PSEUDO_CODON_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter muParam = (Parameter) xo.getElementFirstChild(MU);
            Parameter omegaParam = (Parameter) xo.getElementFirstChild(OMEGA);
            Parameter kappaParam = (Parameter) xo.getElementFirstChild(KAPPA);

            GammaSiteModel siteModel1 = (GammaSiteModel) xo.getElementFirstChild(FIRST_POSITION);
            GammaSiteModel siteModel2 = (GammaSiteModel) xo.getElementFirstChild(SECOND_POSITION);
            GammaSiteModel siteModel3 = (GammaSiteModel) xo.getElementFirstChild(THIRD_POSITION);

            if (!(siteModel1.getSubstitutionModel() instanceof GTR) ||
                    !(siteModel2.getSubstitutionModel() instanceof GTR) ||
                    !(siteModel3.getSubstitutionModel() instanceof GTR)) {
                throw new XMLParseException("Substitution models in " + getParserName() + " elements must be GTRs");
            }

            GTR gtr = (GTR) siteModel1.getSubstitutionModel();
            FrequencyModel freqModel = gtr.getFrequencyModel();

            return new PseudoCodonModel(siteModel1, siteModel2, siteModel3, muParam, omegaParam, kappaParam, freqModel);
        }
        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the Pseudo-Codon model of nucleotide evolution.";
        }

        public Class getReturnType() {
            return PseudoCodonModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(MU, Parameter.class),
                new ElementRule(OMEGA, Parameter.class),
                new ElementRule(KAPPA, Parameter.class),
                new ElementRule(FIRST_POSITION, GammaSiteModel.class),
                new ElementRule(SECOND_POSITION, GammaSiteModel.class),
                new ElementRule(THIRD_POSITION, GammaSiteModel.class)
        };
    };

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the site models
     */
    private GammaSiteModel siteModel1, siteModel2, siteModel3;
    private GTR gtr1, gtr2, gtr3;
    private FrequencyModel frequencyModel;
    private boolean excludeStopCodons;

    protected Parameter muParameter = null;
    protected Parameter omegaParameter = null;
    protected Parameter kappaParameter = null;
}