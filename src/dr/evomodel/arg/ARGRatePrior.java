/*
 * ARGRatePrior.java
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

import dr.evomodel.arg.ARGModel.Node;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

public class ARGRatePrior extends AbstractModelLikelihood {

    public static final String ARG_RATE_PRIOR = "argRatePrior";
    public static final String SIGMA = "sigma";

    private final ARGModel arg;
    private final Parameter logNormalSigma;


    public ARGRatePrior(String name, ARGModel arg, Parameter sigma) {
        super(name);

        this.arg = arg;
        this.logNormalSigma = sigma;

        addModel(arg);
        addVariable(sigma);

    }

    public double[] generateValues() {

        double[] values = new double[arg.getNumberOfPartitions()];

        double sigma = logNormalSigma.getParameterValue(0);

        double oneOverSigma = 1.0 / sigma;

        for (int i = 0; i < values.length; i++) {
            values[i] = MathUtils.nextGamma(oneOverSigma, oneOverSigma);
        }


        return values;
    }

    public double getLogLikelihood() {
        return calculateLogLikelihood();
    }

    public double getAddHastingsRatio(double[] values) {
        return -calculateLogLikelihood(values);
    }


    private double calculateLogLikelihood(double[] values) {
        double logLike = 0;

        double sigma = logNormalSigma.getParameterValue(0);
        double oneOverSigma = 1.0 / sigma;

        for (double d : values) {
            logLike += GammaDistribution.logPdf(d, oneOverSigma, sigma);
        }

        return logLike;
    }

    private double calculateLogLikelihood() {
        double logLike = 0;


        for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
            Node x = (Node) arg.getNode(i);

            if (!x.isRoot() && x.isBifurcation()) {

                double[] values = x.rateParameter.getParameterValues();

                logLike += calculateLogLikelihood(values);
            }
        }


        return logLike;
    }


    public Model getModel() {
        return this;
    }

    public void makeDirty() {

    }

    public String getId() {
        return super.getId();
    }

    public void setId(String id) {
        super.setId(id);
    }

    protected void acceptState() {
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected void handleVariableChangedEvent(Variable variable, int index,
                                               ChangeType type) {
    }

    protected void restoreState() {

    }

    protected void storeState() {
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return null;
        }

        public Class getReturnType() {
            return ARGRatePrior.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return null;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = xo.getAttribute(XMLParser.ID, "");


            Parameter sigma = (Parameter) xo.getChild(Parameter.class);


            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            return new ARGRatePrior(id, arg, sigma);
        }

        public String getParserName() {
            return ARG_RATE_PRIOR;
        }

    };


}
