/*
 * MarginalLikelihoodEstimationGenerator.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.components.marginalLikelihoodEstimation;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.WorkingPriorParsers;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MarginalLikelihoodEstimationGenerator extends BaseComponentGenerator {

    private BeautiOptions beautiOptions = null;

    MarginalLikelihoodEstimationGenerator(final BeautiOptions options) {
        super(options);
        this.beautiOptions = options;
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        if (!component.performMLE && !component.performMLEGSS) {
            return false;
        }

        switch (point) {
            case AFTER_MCMC:
                return true;
        }
        return false;
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);

        switch (point) {
            case AFTER_MCMC:
                writeMLE(writer, component);
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Marginal Likelihood Estimator";
    }

    /**
     * Write the marginalLikelihoodEstimator, pathSamplingAnalysis and steppingStoneSamplingAnalysis blocks.
     *
     * @param writer XMLWriter
     */
    public void writeMLE(XMLWriter writer, MarginalLikelihoodEstimationOptions options) {

        if (options.performMLE) {

            writer.writeComment("Define marginal likelihood estimator (PS/SS) settings");

            List<Attribute> attributes = new ArrayList<Attribute>();
            //attributes.add(new Attribute.Default<String>(XMLParser.ID, "mcmc"));
            attributes.add(new Attribute.Default<Integer>("chainLength", options.mleChainLength));
            attributes.add(new Attribute.Default<Integer>("pathSteps", options.pathSteps));
            attributes.add(new Attribute.Default<String>("pathScheme", options.pathScheme));
            if (!options.pathScheme.equals("linear")) {
                attributes.add(new Attribute.Default<Double>("alpha", options.schemeParameter));
            }

            writer.writeOpenTag("marginalLikelihoodEstimator", attributes);

            writer.writeOpenTag("samplers");
            writer.writeIDref("mcmc", "mcmc");
            writer.writeCloseTag("samplers");

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "pathLikelihood"));
            writer.writeOpenTag("pathLikelihood", attributes);
            writer.writeOpenTag("source");
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, CompoundLikelihoodParser.POSTERIOR);
            writer.writeCloseTag("source");
            writer.writeOpenTag("destination");
            writer.writeIDref(CompoundLikelihoodParser.PRIOR, CompoundLikelihoodParser.PRIOR);
            writer.writeCloseTag("destination");
            writer.writeCloseTag("pathLikelihood");

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "MLELog"));
            attributes.add(new Attribute.Default<Integer>("logEvery", options.mleLogEvery));
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("log", attributes);
            writer.writeIDref("pathLikelihood", "pathLikelihood");
            writer.writeCloseTag("log");

            writer.writeCloseTag("marginalLikelihoodEstimator");

            writer.writeComment("Path sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("pathSamplingAnalysis", attributes);
            writer.writeTag("likelihoodColumn", new Attribute.Default<String>("name", "pathLikelihood.delta"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag("pathSamplingAnalysis");

            writer.writeComment("Stepping-stone sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("steppingStoneSamplingAnalysis", attributes);
            writer.writeTag("likelihoodColumn", new Attribute.Default<String>("name", "pathLikelihood.delta"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag("steppingStoneSamplingAnalysis");

        } else if (options.performMLEGSS) {

            //First define necessary components for the tree working prior
            if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions")) {
                //more general product of exponentials needs to be constructed

                System.err.println("productOfExponentials selected: " + options.choiceTreeWorkingPrior);

                List<Attribute> attributes = new ArrayList<Attribute>();
                attributes.add(new Attribute.Default<String>(XMLParser.ID, "exponentials"));
                attributes.add(new Attribute.Default<String>("fileName", beautiOptions.logFileName));
                attributes.add(new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10));
                attributes.add(new Attribute.Default<String>("parameterColumn", "coalescentEventsStatistic"));
                attributes.add(new Attribute.Default<String>("dimension", "" + (beautiOptions.taxonList.getTaxonCount()-1)));

                writer.writeOpenTag("productOfExponentialsPosteriorMeansLoess", attributes);
                writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.ID, TreeModel.TREE_MODEL), true);
                writer.writeCloseTag("productOfExponentialsPosteriorMeansLoess");

            } else {
                //matching coalescent model has to be constructed
                //getting the coalescent model
                System.err.println("matching coalescent model selected: " + options.choiceTreeWorkingPrior);

                System.err.println(beautiOptions.getPartitionTreePriors().get(0).getNodeHeightPrior());

                //TODO: add the simple parametric coalescent models
                //TODO: if not a simple coalescent model, switch to product of exponentials




            }

            writer.writeComment("Define marginal likelihood estimator (GSS) settings");

            List<Attribute> attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<Integer>("chainLength", options.mleChainLength));
            attributes.add(new Attribute.Default<Integer>("pathSteps", options.pathSteps));
            attributes.add(new Attribute.Default<String>("pathScheme", options.pathScheme));
            if (!options.pathScheme.equals("linear")) {
                attributes.add(new Attribute.Default<Double>("alpha", options.schemeParameter));
            }

            writer.writeOpenTag("marginalLikelihoodEstimator", attributes);

            writer.writeOpenTag("samplers");
            writer.writeIDref("mcmc", "mcmc");
            writer.writeCloseTag("samplers");

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "pathLikelihood"));
            writer.writeOpenTag("pathLikelihood", attributes);
            writer.writeOpenTag("source");
            writer.writeIDref(CompoundLikelihoodParser.POSTERIOR, CompoundLikelihoodParser.POSTERIOR);
            writer.writeCloseTag("source");
            writer.writeOpenTag("destination");
            writer.writeOpenTag("referencePrior");

            ArrayList<Parameter> parameters = beautiOptions.selectParameters();

            for (Parameter param : parameters) {
                System.err.println(param.toString() + "   " + param.priorType.toString());
                //should leave out those parameters set by the coalescent
                if (param.priorType != PriorType.NONE_TREE_PRIOR) {
                    //frequencies is multidimensional, is that automatically dealt with?
                    writer.writeOpenTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR,
                            new Attribute[]{
                                    new Attribute.Default<String>("fileName", beautiOptions.logFileName),
                                    new Attribute.Default<String>("parameterColumn", param.getName()),
                                    new Attribute.Default<String>("burnin", "" + beautiOptions.chainLength*0.10)
                            });
                    writeParameterIdref(writer, param);
                    writer.writeCloseTag(WorkingPriorParsers.NORMAL_REFERENCE_PRIOR);
                }
            }

            if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions")) {
                writer.writeIDref("productOfExponentialsPosteriorMeansLoess", "exponentials");
            } else {
                //TODO: complete this section
            }

            writer.writeCloseTag("referencePrior");
            writer.writeCloseTag("destination");
            writer.writeCloseTag("pathLikelihood");

            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>(XMLParser.ID, "MLELog"));
            attributes.add(new Attribute.Default<Integer>("logEvery", options.mleLogEvery));
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("log", attributes);
            writer.writeIDref("pathLikelihood", "pathLikelihood");
            writer.writeCloseTag("log");

            writer.writeCloseTag("marginalLikelihoodEstimator");

            writer.writeComment("Generalized stepping-stone sampling estimator from collected samples");
            attributes = new ArrayList<Attribute>();
            attributes.add(new Attribute.Default<String>("fileName", options.mleFileName));
            writer.writeOpenTag("generalizedSteppingStoneSamplingAnalysis", attributes);
            writer.writeTag("sourceColumn", new Attribute.Default<String>("name", "pathLikelihood.source"), true);
            writer.writeTag("destinationColumn", new Attribute.Default<String>("name", "pathLikelihood.destination"), true);
            writer.writeTag("thetaColumn", new Attribute.Default<String>("name", "pathLikelihood.theta"), true);
            writer.writeCloseTag("generalizedSteppingStoneSamplingAnalysis");

        }

    }

    private void writeParameterIdref(XMLWriter writer, Parameter parameter) {
        if (parameter.isStatistic) {
            writer.writeIDref("statistic", parameter.getName());
        } else {
            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
        }
    }

}
