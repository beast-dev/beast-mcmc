/*
 * HierarchicalPartitionLikelihood.java
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

import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.xml.*;

public class HierarchicalPartitionLikelihood extends ARGPartitionLikelihood {

    public static final String HIERARCHICAL_PARTITION_LIKELIHOOD = "hierarchicalPartitionLikelihood";

    private Parameter probabilities;

    public HierarchicalPartitionLikelihood(String id, ARGModel arg, Parameter probs) {
        super(id, arg);

        this.probabilities = probs;

        addVariable(probs);
        addModel(arg);
    }

    public double[] generatePartition() {
        double[] partition = new double[getNumberOfPartitionsMinusOne() + 1];

        partition[0] = 0.0;

        for (int i = 0; i < partition.length; i++)
            partition[i] = 0.0;

        while (UniformPartitionLikelihood.arraySum(partition) == 0.0) {
            for (int i = 1; i < partition.length; i++) {
                if (MathUtils.nextDouble() < probabilities.getParameterValue(i - 1)) {
                    partition[i] = 1.0;
                } else {
                    partition[i] = 0.0;
                }
            }
        }


        return partition;
    }

    public double getLogLikelihood(double[] partition) {
        double logLike = 0;

        for (int i = 1; i < partition.length; i++) {
            if (partition[i] == 1.0) {
                logLike += Math.log(probabilities.getParameterValue(i - 1));
            } else {
                logLike += Math.log(1 - probabilities.getParameterValue(i - 1));
            }
        }

//		return 1;

        return logLike;
    }

    protected void acceptState() {
        // nothing to do!
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // i'm lazy

    }

    protected void handleVariableChangedEvent(Variable variable, int index,
                                               ChangeType type) {
        // I'm lazy, so I compute after each step :)
    }

    protected void restoreState() {
        //nothing to restore!
    }

    @Override
    protected void storeState() {
        //nothing to store

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return null;
        }

        public Class getReturnType() {
            return PoissonPartitionLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{

                    new ElementRule(ARGModel.class, false),
            };
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = "";
            if (xo.hasId())
                id = xo.getId();

            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            Parameter values = (Parameter) xo.getChild(Parameter.class);

            if (values.getDimension() != arg.getNumberOfPartitions() - 1) {
                throw new XMLParseException("The dimension of the parameter must equal the number of partitions minus 1 ");
            }

            if (arg.isRecombinationPartitionType()) {
                throw new XMLParseException(ARGModel.TREE_MODEL + " must be of type " + ARGModel.REASSORTMENT_PARTITION);
            }

            return new HierarchicalPartitionLikelihood(id, arg, values);
        }

        public String getParserName() {
            return HIERARCHICAL_PARTITION_LIKELIHOOD;
        }

    };

}
