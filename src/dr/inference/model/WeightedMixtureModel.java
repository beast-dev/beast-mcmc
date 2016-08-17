/*
 * WeightedMixtureModel.java
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
import dr.math.LogTricks;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class WeightedMixtureModel extends AbstractModelLikelihood implements Citable {

    public static final String MIXTURE_MODEL = "mixtureModel";
//    public static final String MIXTURE_WEIGHTS = "weights";
    public static final String NORMALIZE = "normalize";

    public WeightedMixtureModel(List<AbstractModelLikelihood> likelihoodList, Parameter mixtureWeights) {
        super(MIXTURE_MODEL);
        this.likelihoodList = likelihoodList;
        this.mixtureWeights = mixtureWeights;
        for (AbstractModelLikelihood model : likelihoodList) {
            addModel(model);
        }
        addVariable(mixtureWeights);

        StringBuilder sb = new StringBuilder();
        sb.append("Constructing a finite mixture model\n");
        sb.append("\tComponents:\n");
        for (AbstractModelLikelihood model : likelihoodList) {
            sb.append("\t\t\t").append(model.getId()).append("\n");
        }
        sb.append("\tMixing parameter: ").append(mixtureWeights.getId()).append("\n");
        sb.append("\tPlease cite:\n");
        sb.append(Citable.Utils.getCitationString((this)));

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
        double logSum = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < likelihoodList.size(); ++i) {
            double pi = mixtureWeights.getParameterValue(i);
            if (pi > 0.0) {            
                logSum = LogTricks.logSum(logSum,
                        Math.log(pi) + likelihoodList.get(i).getLogLikelihood());
            }
        }   
        return logSum;                
    }

    public void makeDirty() {
    }

    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MIXTURE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter weights = (Parameter) xo.getChild(Parameter.class);
            List<AbstractModelLikelihood> likelihoodList = new ArrayList<AbstractModelLikelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Likelihood)
                    likelihoodList.add((AbstractModelLikelihood) xo.getChild(i));
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

            if (!normalized(weights))
                throw new XMLParseException("Parameter +" + weights.getId() + " must lie on the simplex");

            return new WeightedMixtureModel(likelihoodList, weights);
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
                new ElementRule(Likelihood.class,2,Integer.MAX_VALUE),
                new ElementRule(Parameter.class)
        };
    };


    private final Parameter mixtureWeights;
    List<AbstractModelLikelihood> likelihoodList;


    public static void main(String[] args) {

        final double l1 = -10;
        final double l2 = -2;

        AbstractModelLikelihood like1 = new AbstractModelLikelihood("dummy") {


            public Model getModel() {
                return null;
            }

            public double getLogLikelihood() {
                return l1;
            }

            public void makeDirty() {
            }

            public String prettyName() {
                return null;
            }

            public boolean isUsed() {
                return false;
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) {
            }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
            }

            @Override
            protected void storeState() {
            }

            @Override
            protected void restoreState() {
            }

            @Override
            protected void acceptState() {
            }

            public void setUsed() {
            }

            public LogColumn[] getColumns() {
                return new LogColumn[0];
            }

            public String getId() {
                return null;
            }

            public void setId(String id) {
            }
        };

        AbstractModelLikelihood like2 = new AbstractModelLikelihood("dummy") {

            public Model getModel() {
                return null;
            }

            public double getLogLikelihood() {
                return l2;
            }

            public void makeDirty() {
            }

            public String prettyName() {
                return null;
            }

            public boolean isUsed() {
                return false;
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) {
            }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
            }

            @Override
            protected void storeState() {
            }

            @Override
            protected void restoreState() {
            }

            @Override
            protected void acceptState() {
            }

            public void setUsed() {
            }

            public LogColumn[] getColumns() {
                return new LogColumn[0];
            }

            public String getId() {
                return null;
            }

            public void setId(String id) {                
            }
        };

        List<AbstractModelLikelihood> likelihoodList = new ArrayList<AbstractModelLikelihood>();
        likelihoodList.add(like1);
        likelihoodList.add(like2);

        Parameter weights = new Parameter.Default(2);
        double p1 = 0.05;
        weights.setParameterValue(0, p1);
        weights.setParameterValue(1, 1.0 - p1);

        WeightedMixtureModel mixture = new WeightedMixtureModel(likelihoodList, weights);
        System.err.println("getLogLikelihood() = " + mixture.getLogLikelihood());

        double test = Math.log(p1 * Math.exp(l1) + (1.0 - p1) * Math.exp(l2));
        System.err.println("correct            = " + test);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MISC;
    }

    @Override
    public String getDescription() {
        return "Weighted mixture model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_MIXTURE_2012);
    }
}
