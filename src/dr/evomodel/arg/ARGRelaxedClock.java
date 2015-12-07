/*
 * ARGRelaxedClock.java
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

package dr.evomodel.arg;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.arg.ARGModel.Node;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.xml.*;

public class ARGRelaxedClock extends AbstractBranchRateModel {

    public static final String ARG_LOCAL_CLOCK = "argLocalClock";
    public static final String PARTITION = "partition";

    private Parameter globalRateParameter;

    private ARGModel arg;
    private int partition;


    public ARGRelaxedClock(String name) {
        super(name);
    }

    public ARGRelaxedClock(String name, ARGModel arg, int partition, Parameter rate) {
        super(name);

        this.arg = arg;
        this.partition = partition;

        globalRateParameter = rate;

        addModel(arg);
        addVariable(rate);
    }

    protected void acceptState() {

    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        //do nothing
    }


    protected void handleVariableChangedEvent(Variable variable, int index, ChangeType type) {
        //do nothing
    }


    protected void restoreState() {

    }


    protected void storeState() {

    }

    public double getBranchRate(Tree tree, NodeRef nodeRef) {

        Node treeNode = (Node) nodeRef;
        Node argNode = (Node) treeNode.mirrorNode;


        return globalRateParameter.getParameterValue(0) * argNode.getRate(partition);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return null;
        }

        public Class getReturnType() {

            return ARGRelaxedClock.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {

            return null;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            int partition = xo.getAttribute(PARTITION, 0);

            Parameter rate = (Parameter) xo.getChild(Parameter.class);

            return new ARGRelaxedClock("", arg, partition, rate);
        }

        public String getParserName() {
            return ARG_LOCAL_CLOCK;
        }

    };

}
