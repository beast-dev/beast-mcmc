/*
 * HierarchicalModelComponentGenerator.java
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

package dr.app.beauti.components.hpm;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.util.XMLWriter;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.ParameterParser;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.NormalGammaPrecisionGibbsOperator;
import dr.inference.operators.NormalNormalMeanGibbsOperator;
import dr.inferencexml.distribution.DistributionLikelihoodParser;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.inferencexml.distribution.NormalDistributionModelParser;
import dr.inferencexml.distribution.PriorParsers;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @version $Id$
 */
public class HierarchicalModelComponentGenerator extends BaseComponentGenerator {

    public static final String MEAN_SUFFIX = ".mean";
    public static final String PRECISION_SUFFIX = ".precision";
    public static final String MODEL_SUFFIX = ".model";
    public static final String HPM_SUFFIX = ".hpm";
    public static final String MEAN_PRIOR_SUFFIX = ".prior.mean";
    public static final String PRECISION_PRIOR_SUFFIX = ".prior.precision";

    public HierarchicalModelComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);

        switch (point) {
            case AFTER_TREE_LIKELIHOOD:
            case IN_OPERATORS:
            case IN_MCMC_PRIOR:
            case IN_FILE_LOG_PARAMETERS:
                return !comp.isEmpty();
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        HierarchicalModelComponentOptions comp = (HierarchicalModelComponentOptions)
                options.getComponentOptions(HierarchicalModelComponentOptions.class);

        switch (point) {
            case AFTER_TREE_LIKELIHOOD:
                generateDistributions(comp.getHPMList(), writer);
                break;
            case IN_OPERATORS:
                generateOperators(comp.getHPMList(), writer);
                break;
            case IN_MCMC_PRIOR:
                generatePriors(comp.getHPMList(), writer);
                break;
            case IN_FILE_LOG_PARAMETERS:
                generateLogs(comp.getHPMList(), writer);
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }
    }

    private void generateOperators(List<HierarchicalPhylogeneticModel> hpmList, XMLWriter writer) {
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            generateOperators(hpm, writer);
        }
    }

    private void generateOperators(HierarchicalPhylogeneticModel hpm, XMLWriter writer) {

        // Generate Normal-normal operator on mean
        writer.writeOpenTag(NormalNormalMeanGibbsOperator.OPERATOR_NAME, getOperatorAttributes());

        writer.writeOpenTag(NormalNormalMeanGibbsOperator.LIKELIHOOD);
        writer.writeIDref(DistributionLikelihoodParser.DISTRIBUTION, getDistributionName(hpm));
        writer.writeCloseTag(NormalNormalMeanGibbsOperator.LIKELIHOOD);

        writer.writeOpenTag(NormalNormalMeanGibbsOperator.PRIOR);
        writer.writeIDref(PriorParsers.NORMAL_PRIOR, getMeanPriorName(hpm));
        writer.writeCloseTag(NormalNormalMeanGibbsOperator.PRIOR);

        writer.writeCloseTag(NormalNormalMeanGibbsOperator.OPERATOR_NAME);

        // Generate Gamma-normal operator on precision
        writer.writeOpenTag(NormalGammaPrecisionGibbsOperator.OPERATOR_NAME, getOperatorAttributes());

        writer.writeOpenTag(NormalGammaPrecisionGibbsOperator.LIKELIHOOD);
        writer.writeIDref(DistributionLikelihoodParser.DISTRIBUTION, getDistributionName(hpm));
        writer.writeCloseTag(NormalGammaPrecisionGibbsOperator.LIKELIHOOD);

        writer.writeOpenTag(NormalGammaPrecisionGibbsOperator.PRIOR);
        writer.writeIDref(PriorParsers.GAMMA_PRIOR, getPrecisionPriorName(hpm));
        writer.writeCloseTag(NormalGammaPrecisionGibbsOperator.PRIOR);

        writer.writeCloseTag(NormalGammaPrecisionGibbsOperator.OPERATOR_NAME);
    }

    private void generateLogs(List<HierarchicalPhylogeneticModel> hpmList, XMLWriter writer) {
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            generateLog(hpm, writer);
        }
    }

    private void generateLog(HierarchicalPhylogeneticModel hpm, XMLWriter writer) {
        for (Parameter parameter : hpm.getConditionalParameterList()) {
            writer.writeIDref(ParameterParser.PARAMETER, parameter.getName());
        }
    }

    private void generatePriors(List<HierarchicalPhylogeneticModel> hpmList, XMLWriter writer) {
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            generatePrior(hpm, writer);
        }
    }

    private void generatePrior(HierarchicalPhylogeneticModel hpm, XMLWriter writer) {
        writer.writeIDref(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD, getDistributionName(hpm));
        writer.writeIDref(PriorParsers.NORMAL_PRIOR, getMeanPriorName(hpm));
        writer.writeIDref(PriorParsers.GAMMA_PRIOR, getPrecisionPriorName(hpm));
    }

    private void generateDistributions(List<HierarchicalPhylogeneticModel> hpmList, XMLWriter writer) {
        for (HierarchicalPhylogeneticModel hpm : hpmList) {
            generateDistribution(hpm, writer);
            generateNormalAndGammaPrior(hpm, writer);
        }
    }

    private String getModelTagName(HierarchicalPhylogeneticModel hpm) {
        switch (hpm.getPriorType()) {
            case NORMAL_HPM_PRIOR:
                return NormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL;
            case LOGNORMAL_HPM_PRIOR:
                return LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL;
            default:
        }
        throw new RuntimeException("Unimplemented HPM prior type");
    }

    private String getDistributionName(HierarchicalPhylogeneticModel hpm) {
        return hpm.getName() + HPM_SUFFIX;
    }

    private String getModelName(HierarchicalPhylogeneticModel hpm) {
        return hpm.getName() + MODEL_SUFFIX;
    }

    private Attribute[] getOperatorAttributes() {
        return new Attribute[] {
                new Attribute.Default<Double>(MCMCOperator.WEIGHT, 1.0),
        };
    }

    private Attribute[] getModelAttributes(HierarchicalPhylogeneticModel hpm) {
        switch (hpm.getPriorType()){
            case NORMAL_HPM_PRIOR:
                return new Attribute[] {
                        new Attribute.Default<String>(XMLParser.ID, getModelName(hpm)),
                };
            case LOGNORMAL_HPM_PRIOR:
                return new Attribute[] {
                        new Attribute.Default<String>(XMLParser.ID, getModelName(hpm)),
                        //new Attribute.Default<Boolean>(LogNormalDistributionModelParser.MEAN_IN_REAL_SPACE, false),
                };
            default:
        }
        throw new RuntimeException("Unimplemented HPM prior type");
    }

    private Attribute[] getMeanPriorAttributes(HierarchicalPhylogeneticModel hpm) {
        return new Attribute[] {
                new Attribute.Default<String>(XMLParser.ID, getMeanPriorName(hpm)),
                new Attribute.Default<Double>(PriorParsers.MEAN, hpm.getConditionalParameterList().get(0).mean),
                new Attribute.Default<Double>(PriorParsers.STDEV, hpm.getConditionalParameterList().get(0).stdev),
        };
    }

    private String getMeanPriorName(HierarchicalPhylogeneticModel hpm) {
        return hpm.getName() + MEAN_PRIOR_SUFFIX;
    }

    private String getPrecisionPriorName(HierarchicalPhylogeneticModel hpm) {
        return hpm.getName() + PRECISION_PRIOR_SUFFIX;
    }

    private Attribute[] getPrecisionPriorAttributes(HierarchicalPhylogeneticModel hpm) {
        return new Attribute[] {
                new Attribute.Default<String>(XMLParser.ID, getPrecisionPriorName(hpm)),
                new Attribute.Default<Double>(PriorParsers.SHAPE, hpm.getConditionalParameterList().get(1).shape),
                new Attribute.Default<Double>(PriorParsers.SCALE, hpm.getConditionalParameterList().get(1).scale),
                new Attribute.Default<Double>(PriorParsers.OFFSET, 0.0),
        };
    }

    private void generateDistribution(HierarchicalPhylogeneticModel hpm, XMLWriter writer) {
      
        writer.writeOpenTag(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD,
                new Attribute[] {
                        new Attribute.Default<String>(XMLParser.ID, getDistributionName(hpm)),
        });

        // Add parameters as data
        writer.writeOpenTag(DistributionLikelihoodParser.DATA);
        for (Parameter parameter : hpm.getArgumentParameterList())  {
            writeParameterRef(parameter.getName(), writer);
        }
        writer.writeCloseTag(DistributionLikelihoodParser.DATA);

        // Add HPM model
        writer.writeOpenTag(DistributionLikelihoodParser.DISTRIBUTION);
        writer.writeOpenTag(getModelTagName(hpm), getModelAttributes(hpm));
        
        writeParameter(NormalDistributionModelParser.MU,
                hpm.getConditionalParameterList().get(0).getName(), 1, hpm.getConditionalParameterList().get(0).getInitial(),
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, writer);
        writeParameter(NormalDistributionModelParser.PREC,
                hpm.getConditionalParameterList().get(1).getName(), 1, hpm.getConditionalParameterList().get(1).getInitial(),
                0.0, Double.POSITIVE_INFINITY, writer);
        writer.writeCloseTag(getModelTagName(hpm));

        writer.writeCloseTag(DistributionLikelihoodParser.DISTRIBUTION);

        writer.writeCloseTag(DistributionLikelihood.DISTRIBUTION_LIKELIHOOD);
    }

    private void generateNormalAndGammaPrior(HierarchicalPhylogeneticModel hpm, XMLWriter writer) {
        // Normal prior on mean
        writer.writeOpenTag(PriorParsers.NORMAL_PRIOR, getMeanPriorAttributes(hpm));
        writer.writeIDref(ParameterParser.PARAMETER, hpm.getConditionalParameterList().get(0).getName());
        writer.writeCloseTag(PriorParsers.NORMAL_PRIOR);

        // Gamma prior on precision
        writer.writeOpenTag(PriorParsers.GAMMA_PRIOR, getPrecisionPriorAttributes(hpm));
        writer.writeIDref(ParameterParser.PARAMETER, hpm.getConditionalParameterList().get(1).getName());
        writer.writeCloseTag(PriorParsers.GAMMA_PRIOR);                
    }

    protected String getCommentLabel() {
        return "Hierarchical phylogenetic models";
    }

}