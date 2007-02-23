/*
 * YuleModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;
//import dr.evomodel.tree.ARGModel;


/**
 * This class contains methods that describe a Yule speciation with Poisson HGT event model.
 *
 * @author Marc Suchard
 */

public class MrBayesDefaultModel extends SpeciationModel {

    public static final String MRBAYES_DEFAULT_MODEL = "mrbayesDefaultModel";
    public static String BIRTH_RATE = "birthRate";
    public static final String HGT_RATE = "hgtRate";


    public MrBayesDefaultModel(Parameter birthRateParameter, Parameter hgtRateParameter, Type units) {

        super(MRBAYES_DEFAULT_MODEL, units);

        this.birthRateParameter = birthRateParameter;
        this.hgtRateParameter = hgtRateParameter;
        addParameter(birthRateParameter);
        addParameter(hgtRateParameter);
        birthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        hgtRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

    }

    public double getBirthRate() {
        return birthRateParameter.getParameterValue(0);
    }

    public double getHgtRate() {
        return hgtRateParameter.getParameterValue(0);
    }

    public void setBirthRate(double birthRate) {

        birthRateParameter.setParameterValue(0, birthRate);
    }

    //
    // functions that define a speciation model
    //
    public double logTreeProbability(int taxonCount) {
        if (true)
            throw new RuntimeException("Why was 'logTreeProbability' called?");
        return 0.0;
    }

    public double logReassortmentProbability(Tree arg) {
        return 0.0;
    }

    //
    // functions that define a speciation model
    //
    public double logNodeProbability(Tree tree, NodeRef node) {
        if (tree.getRoot() != node)
            //  tree.set
            return 0.0;

//        ARGModel arg = (ARGModel)tree;

//        if (arg.isReassortment(node))
//            return 0.0;

//        double nodeHeight = tree.getNodeHeight(node);
        double rootHeight = tree.getNodeHeight(node);

        double lambda = getBirthRate();
//        System.err.println(getBirthRate());
        //return Math.log((lambda * Math.exp(-lambda * nodeHeight)) / (1 - Math.exp(-lambda * rootHeight)));
        return -lambda * rootHeight;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    /**
     * Parses an element from an DOM document into a SpeciationModel. Recognises
     * YuleModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MRBAYES_DEFAULT_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            //int units = XMLParser.Utils.getUnitsAttr(xo);
            Type units = XMLParser.Utils.getUnitsAttr(xo);

            XMLObject cxo = (XMLObject) xo.getChild(BIRTH_RATE);
            Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);
            cxo = (XMLObject) xo.getChild(HGT_RATE);
            Parameter hgtParameter = (Parameter) cxo.getChild(Parameter.class);

            return new MrBayesDefaultModel(brParameter, hgtParameter, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model of a Yule with Horizontal Gene Transfer process.";
        }

        public Class getReturnType() {
            return YuleModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(BIRTH_RATE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class), new ElementRule(Parameter.class)}),
                XMLUnits.SYNTAX_RULES[0]
        };
    };


    //Protected stuff
    private Parameter birthRateParameter;
    private Parameter hgtRateParameter;
}
