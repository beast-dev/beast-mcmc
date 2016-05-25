/*
 * PoissonPartitionLikelihood.java
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
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.inferencexml.distribution.PriorParsers;
import dr.math.MathUtils;
import dr.math.Poisson;
import dr.math.distributions.PoissonDistribution;
import dr.xml.*;
import jebl.math.Binomial;

public class PoissonPartitionLikelihood extends ARGPartitionLikelihood {

    public static final String POISSON_PARTITION_LIKELIHOOD = "poissonPartitionLikelihood";

    PoissonDistribution pd;
    double mean;

    public PoissonPartitionLikelihood(String id, ARGModel arg, double mean) {
        super(id, arg);

        pd = new PoissonDistribution(mean);
        this.mean = mean;
    }

    public double getLogLikelihood(double[] partition) {
        if ((getNumberOfPartitionsMinusOne() + 1) % 2 == 0) {
            return getEvenLogLikelihood(partition);
        }

        return getOddLogLikelihood(partition);

    }

    private double getEvenLogLikelihood(double[] partition) {
        int numberOfZeros = 0;
        int numberOfOnes = 0;

        for (double d : partition) {
            assert d == 0.0 || d == 1.0;

            if (d == 0.0)
                numberOfZeros++;
            else
                numberOfOnes++;
        }

        double poissonValue = (double) Math.min(numberOfZeros, numberOfOnes);

        double logLike = pd.logPdf(poissonValue);

        if (poissonValue < partition.length / 2) {
            return logLike - Math.log(Binomial.choose(partition.length, poissonValue));
        } else {
            return logLike - Math.log(Binomial.choose(partition.length - 1, poissonValue));
        }
    }

    private double getOddLogLikelihood(double[] partition) {
        int numberOfZeros = 0;
        int numberOfOnes = 0;

        for (double d : partition) {
            assert d == 0.0 || d == 1.0;

            if (d == 0.0)
                numberOfZeros++;
            else
                numberOfOnes++;

        }

        double poissonValue = (double) Math.min(numberOfZeros, numberOfOnes);

        return pd.logPdf(poissonValue) -
                Math.log(Binomial.choose(partition.length, poissonValue));
    }

    public double[] generatePartition() {
        int lengthDividedByTwo = (getNumberOfPartitionsMinusOne() + 1) / 2;

        int value = 0;

        while (value < 1 || value > lengthDividedByTwo) {
            value = Poisson.nextPoisson(mean);
        }

        int[] x = new int[getNumberOfPartitionsMinusOne() + 1];

        for (int i = 0; i < value; i++) {
            x[i] = 1;
        }

        MathUtils.permute(x);

        if (x[0] == 1) {
            for (int i = 0; i < x.length; i++) {
                if (x[i] == 1) {
                    x[i] = 0;
                } else {
                    x[i] = 1;
                }
            }
        }

        double[] rValue = new double[x.length];

        for (int i = 0; i < rValue.length; i++) {
            rValue[i] = x[i];
        }

        return rValue;
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
                    AttributeRule.newDoubleRule(PriorParsers.MEAN, false),
                    new ElementRule(ARGModel.class, false),
            };
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = "";
            if (xo.hasId())
                id = xo.getId();

            double mean = xo.getDoubleAttribute(PriorParsers.MEAN);
            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            if (arg.isRecombinationPartitionType()) {
                throw new XMLParseException(ARGModel.TREE_MODEL + " must be of type " + ARGModel.REASSORTMENT_PARTITION);
            }

            return new PoissonPartitionLikelihood(id, arg, mean);
        }

        public String getParserName() {
            return POISSON_PARTITION_LIKELIHOOD;
        }

    };

    protected void acceptState() {
        //nothing to do
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        //has no submodels
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                               ChangeType type) {
        //has no parameters
    }

    @Override
    protected void restoreState() {
        //nothing to restore
    }

    @Override
    protected void storeState() {
        // nothing to store
    }


}
