/*
 * MixtureModelLikelihood.java
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

package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 * @author Alexander V. Alekseyenko
 */
public class MixtureModelLikelihood extends AbstractModelLikelihood implements Citable, Loggable {

    public static final String MIXTURE_MODEL_ALIAS = "integratedMixtureModel";
    public static final String MIXTURE_MODEL = "mixtureModelLikelihood";
    //    public static final String MIXTURE_WEIGHTS = "weights";
    public static final String NORMALIZE = "normalize";

    public MixtureModelLikelihood(List<Likelihood> likelihoodList, Parameter weights) {
        super(MIXTURE_MODEL);
        this.likelihoodList = likelihoodList;
        this.mixtureWeights = weights;
        for (Likelihood model : likelihoodList) {
            if (model.getModel() != null) {
                addModel(model.getModel());
            }
        }
        addVariable(mixtureWeights);

        StringBuilder sb = new StringBuilder();
        sb.append("Constructing a finite mixture model\n");
        sb.append("\tComponents:\n");
        for (Likelihood like : likelihoodList) {
            Model model = like.getModel();
            sb.append("\t\t\t").append(
                    model != null ?
                    like.getModel().getId() : "anonymous"
            ).append("\n");
        }
//        sb.append("\tMixing parameter: ").append(mixtureWeights.getId()).append("\n");
        sb.append("\tPlease cite:\n");
        sb.append(Utils.getCitationString((this)));

        Logger.getLogger("dr.inference.model").info(sb.toString());
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return getLogLikelihoodSum();
    }

    private double getLogLikelihoodSum() {
        double logSum = Double.NEGATIVE_INFINITY;
        double bad = 0;

        double[] weights;
        if(useParameter)
            weights = MathUtils.getNormalized(mixtureWeights.getParameterValues());
        else{
            weights = new double[likelihoodList.size()];
            for(int i=0; i<likelihoodList.size(); ++i) weights[i]=1.0/likelihoodList.size();
        }

        for (int i = 0; i < likelihoodList.size(); ++i) {
            double pi = weights[i];
            if (pi > 0.0) {
                logSum = LogTricks.logSum(logSum,
                        Math.log(pi) + likelihoodList.get(i).getLogLikelihood());
            }
            bad += likelihoodList.get(i).getLogLikelihood() * pi;
        }
        if (powerPrior) {
            return bad;
        } else {
            return logSum;
        }
    }

    private double getWeight(final int dim) {
        if (useParameter) {
            return mixtureWeights.getParameterValue(dim);
        } else {
            return 1.0 / likelihoodList.size();
        }
    }

    public void makeDirty() {
        // Do nothing
    }

    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[likelihoodList.size()];
        for (int i = 0; i < likelihoodList.size(); ++i) {
            columns[i] = new MixtureColumn(MIXTURE_MODEL, i);
        }
        return columns;
    }

    private class MixtureColumn extends NumberColumn {

        public MixtureColumn(String label, int dim) {
            super(label);
            this.dim = dim;
        }

        @Override
        public double getDoubleValue() {
            double logSum = getLogLikelihoodSum();
            double logLike = likelihoodList.get(dim).getLogLikelihood() +  Math.log(getWeight(dim));



            double x =  logLike - logSum;

            if (inProbSpace) {
                x = Math.exp(x);
            }

//            System.err.println(logLike + " : " + logSum + " " + dim + " " + x);
//            System.exit(-1);
            return x;
        }

        private final int dim;
        private final boolean inProbSpace = true;
    }


    private static final boolean useParameter = true;
    private static final boolean powerPrior = false;

    public static XMLObjectParser PARSER_ALIAS = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            try{
                return ((AbstractXMLObjectParser)PARSER).parseXMLObject(xo);
            }
            catch(XMLParseException e){
                throw(e);
            }
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return PARSER.getSyntaxRules();
        }


        public String getParserDescription() {
            return PARSER.getParserDescription();
        }


        public Class getReturnType() {
            return PARSER.getReturnType();
        }

        public String getParserName() {
            return MIXTURE_MODEL_ALIAS;
        }
    };

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MIXTURE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter weights = (Parameter) xo.getChild(Parameter.class);
            List<Likelihood> likelihoodList = new ArrayList<Likelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Likelihood)
                    likelihoodList.add((Likelihood) xo.getChild(i));
            }

            if (weights.getDimension() != likelihoodList.size()) {
                throw new XMLParseException("Dim of " + weights.getId() + " does not match the number of likelihoods");
            }

            if (xo.hasAttribute(NORMALIZE)) {
                if (xo.getBooleanAttribute(NORMALIZE)) {
                    double sum = 0;
                    for (int i = 0; i < weights.getDimension(); i++)
                        sum += weights.getParameterValue(i);
                    for (int i = 0; i < weights.getDimension(); i++)
                        weights.setParameterValue(i, weights.getParameterValue(i) / sum);
                }
            }

            return new MixtureModelLikelihood(likelihoodList, weights);
        }

        private boolean normalized(Parameter p) {
            double sum = 0;
            for (int i = 0; i < p.getDimension(); i++)
                sum += p.getParameterValue(i);
            return (sum == 1.0);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a finite mixture of likelihood models.";
        }

        public Class getReturnType() {
            return CompoundModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(NORMALIZE, true),
                new ElementRule(Likelihood.class, 2, Integer.MAX_VALUE),
                new ElementRule(Parameter.class)
        };
    };
    private final Parameter mixtureWeights;
    List<Likelihood> likelihoodList;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MISC;
    }

    @Override
    public String getDescription() {
        return "Mixture model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_MIXTURE_2012);
    }

}
