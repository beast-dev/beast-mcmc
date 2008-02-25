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


/**
 * The Yule model.
 * <p/>
 * The trouble seems to be that Nee (2001) and Yang & Rannala (1997) are both
 * avoiding assuming a prior on t_1 (which becomes a hyperprior in the context
 * of Bayesian phylogenetics). That is, what do you know about t_1 before
 * you've seen the sample, and don't even know the number of species? It seems
 * an improper uniform prior is the obvious choice, and that is essentially
 * what you need to go from an expression like Nee 2001, eq (3) to a
 * probability for t_1. Simplest possible case: pure birth process with rate v,
 * and two tips.
 * <p/>
 * P(tree with 2 tips|root at time t) = exp(-2vt)
 * <p/>
 * This is the likelihood
 * <p/>
 * L(root at time t|tree with 2 tips)
 * <p/>
 * but to be a probability density, it needs normalising, so with a uniform
 * prior on t, this is
 * <p/>
 * 2v exp(-2vt).
 * <p/>
 * -- comments by Graham Jones, pers. comms.
 *
 * @author Alexei Drummond
 * @author Roald Forsberg
 */
public class YuleModel extends SpeciationModel {

    public static final String YULE_MODEL = "yuleModel";
    public static String BIRTH_RATE = "birthRate";


    public YuleModel(Parameter birthRateParameter, Type units) {

        super(YULE_MODEL, units);

        this.birthRateParameter = birthRateParameter;
        addParameter(birthRateParameter);
        birthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }

    public double getBirthRate() {
        return birthRateParameter.getParameterValue(0);
    }

    public void setBirthRate(double birthRate) {

        birthRateParameter.setParameterValue(0, birthRate);
    }

    //
    // functions that define a speciation model
    //
    public double logTreeProbability(int taxonCount) {
        return 0.0;
    }

    //
    // functions that define a speciation model
    //
    public double logNodeProbability(Tree tree, NodeRef node) {

        double nodeHeight = tree.getNodeHeight(node);
        final double lambda = getBirthRate();

        double logP = 0;


        if (tree.isRoot(node)) {
            // see Appendix 1 of Nee (2001) paper for discussion about why we double this
            // nodeHeight for the root.
            nodeHeight *= 2;
        } else {
            // see Appendix 1 of Nee (2001) paper for discussion about why we leave off
            // this contribution for the last internode
            logP += Math.log(lambda);
        }
        logP += -lambda * nodeHeight;

        return logP;
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return false;
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
            return YULE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLParser.Utils.getUnitsAttr(xo);

            XMLObject cxo = (XMLObject) xo.getChild(BIRTH_RATE);
            Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);

            return new YuleModel(brParameter, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A speciation model of a Yule process.";
        }

        public Class getReturnType() {
            return YuleModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(BIRTH_RATE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                XMLUnits.SYNTAX_RULES[0]
        };
    };


    //Protected stuff
    Parameter birthRateParameter;

}
