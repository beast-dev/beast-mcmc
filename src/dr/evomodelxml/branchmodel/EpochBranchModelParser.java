/*
 * EpochBranchModelParser.java
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

package dr.evomodelxml.branchmodel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.branchratemodel.RateEpochBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 */
public class EpochBranchModelParser extends AbstractXMLObjectParser {

    public static final String EPOCH_BRANCH_MODEL = "epochBranchModel";
    public static final String EPOCH = "epoch";
    public static final String TRANSITION_TIME = "transitionTime";

    public String getParserName() {
        return EPOCH_BRANCH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("\nUsing multi-epoch branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel ancestralSubstitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        List<Epoch> epochs = new ArrayList<Epoch>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {
                XMLObject xoc = (XMLObject) xo.getChild(i);
                if (xoc.getName().equals(EPOCH)) {

                    Parameter tt = null;

                    if (xoc.hasAttribute(TRANSITION_TIME)) {
                        double t = xoc.getAttribute(TRANSITION_TIME, 0.0);
                        tt = new Parameter.Default(1, t);
                    }

                    SubstitutionModel s = (SubstitutionModel) xoc.getChild(SubstitutionModel.class);

                    if (xoc.hasChildNamed(TRANSITION_TIME)) {
                        if (tt != null) {
                            throw new XMLParseException("An epoch cannot have a transitionTime attribute and a parameter");
                        }

                        tt = (Parameter) xoc.getElementFirstChild(TRANSITION_TIME);
                    }
                    epochs.add(new Epoch(s, tt));
                }
            }
        }

        Collections.sort(epochs);
        List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
        CompoundParameter transitionTimes = new CompoundParameter("epochTimes");

        for (Epoch epoch : epochs) {
            substitutionModels.add(epoch.substitutionModel);
            transitionTimes.addParameter(epoch.timeParameter);
        }

        substitutionModels.add(ancestralSubstitutionModel);

        return new EpochBranchModel(tree, substitutionModels, transitionTimes);
    }

    class Epoch implements Comparable<Object> {

        private final double transitionTime;
        private final SubstitutionModel substitutionModel;
        private final Parameter timeParameter;

        public Epoch(SubstitutionModel substitutionModel, Parameter timeParameter) {
            this.transitionTime = timeParameter.getParameterValue(0);
            this.substitutionModel = substitutionModel;
            this.timeParameter = timeParameter;
        }

        public int compareTo(Object o) {
            return Double.compare(transitionTime, ((Epoch) o).transitionTime);
        }

    }
    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a branch model which has multiple epoch. " +
                        "All branches (or portions of them) have the same substitution " +
                        "model within a given epoch. If parameters are used to sample " +
                        "transition times, these must be kept in ascending order by judicious " +
                        "use of bounds or priors.";
    }

    public Class<RateEpochBranchRateModel> getReturnType() {
        return RateEpochBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class, "The tree across which the epochs will be assigned"),
            new ElementRule(SubstitutionModel.class, "The ancestral substitution model after the last epoch"),
            new ElementRule(EPOCH,
                    new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(TRANSITION_TIME, true, "The time of transition between this epoch and the previous one"),
                            new ElementRule(SubstitutionModel.class, "The substitution model for this epoch"),
                            new ElementRule(TRANSITION_TIME, Parameter.class, "The transition time parameter for this epoch", true)
                    }, "An epoch that lasts until transitionTime",
                    1, Integer.MAX_VALUE
            ),
    };

}
