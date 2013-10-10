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

import dr.app.beauti.components.sequenceerror.SequenceErrorModelComponentOptions;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.types.SequenceErrorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.HypermutantAlignment;
import dr.evomodel.treelikelihood.HypermutantErrorModel;
import dr.evomodelxml.treelikelihood.SequenceErrorModelParser;
import dr.evoxml.HypermutantAlignmentParser;
import dr.inference.model.ParameterParser;
import dr.inference.model.StatisticParser;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MarginalLikelihoodEstimationGenerator extends BaseComponentGenerator {

    MarginalLikelihoodEstimationGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        MarginalLikelihoodEstimationOptions component = (MarginalLikelihoodEstimationOptions) options.getComponentOptions(MarginalLikelihoodEstimationOptions.class);


        if (!component.performMLE) {
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

            writer.writeComment("Define marginal likelihood estimator settings");

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

        }

    }

}
