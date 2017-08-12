/*
 * DolloComponentGenerator.java
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

package dr.app.beauti.components.dollo;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodelxml.MSSD.ALSSiteModelParser;
import dr.oldevomodelxml.MSSD.ALSTreeLikelihoodParser;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.oldevomodelxml.sitemodel.GammaSiteModelParser;
import dr.oldevomodelxml.substmodel.MutationDeathModelParser;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.MutationDeathTypeParser;
import dr.evoxml.SitePatternsParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Marc Suchard
 * @version $Id$
 */

public class DolloComponentGenerator extends BaseComponentGenerator {

	protected DolloComponentGenerator(BeautiOptions options) {
		super(options);
	}

	public boolean usesInsertionPoint(InsertionPoint point) {

		DolloComponentOptions component = (DolloComponentOptions) options
				.getComponentOptions(DolloComponentOptions.class);

        if (!component.isActive()) {
            return false;
        }

		switch (point) {
        case AFTER_TAXA:
        case AFTER_SUBSTITUTION_MODEL:
        case IN_FILE_LOG_PARAMETERS:
//        case IN_MCMC_LIKELIHOOD:
//        case AFTER_TREE_LIKELIHOOD:
//		case IN_OPERATORS:
//		case IN_TREES_LOG:
//		case AFTER_TREES_LOG:
//		case AFTER_MCMC:
			return true;
		default:
			return false;
		}

	}// END: usesInsertionPoint

	protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {

		DolloComponentOptions component = (DolloComponentOptions) options
				.getComponentOptions(DolloComponentOptions.class);

        if (!includeStochasticDollo()) {
            return;
        }

		switch (point) {
        case AFTER_TAXA:
            writeDataType(writer);
            break;
        case AFTER_SUBSTITUTION_MODEL:
            writeDolloSubstitutionModels(writer, component);
            break;
//        case AFTER_TREE_LIKELIHOOD:
//            writeDolloTreeLikelihoods(writer, component);
//            break;
		case IN_OPERATORS:
			break;
//        case IN_MCMC_LIKELIHOOD:
//            writeDolloTreeLikelihoodReferences(writer);
//            break;
            
        case IN_SCREEN_LOG:
            writeScreenLogEntries(writer, component);
            break;            
		case IN_FILE_LOG_PARAMETERS:
            writeLog(writer, component);
			break;
		case IN_TREES_LOG:
			break;
		case AFTER_TREES_LOG:
			break;
		case AFTER_MCMC:
			break;
		default:
			throw new IllegalArgumentException(
					"This insertion point is not implemented for "
							+ this.getClass().getName());
		}

	}// END: generate

    private void writeScreenLogEntries(XMLWriter writer, DolloComponentOptions component) {
    }

//    private void writeDolloTreeLikelihoodReferences(XMLWriter writer) {
//
//        for (AbstractPartitionData partition : options.dataPartitions) {
//            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
//            String prefix = partition.getName();
//            if (model.isDolloModel()) {
//                writer.writeIDref(ALSTreeLikelihoodParser.LIKE_NAME,
//                        prefix + ALSTreeLikelihoodParser.LIKE_NAME);
//            }
//        }
//    }

    private void writeLog(XMLWriter writer, DolloComponentOptions component) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model.isDolloModel()) {
                String prefix = partition.getName() + ".";
                writer.writeIDref(ParameterParser.PARAMETER,
                    component.getOptions().getParameter(prefix + DolloComponentOptions.DEATH_RATE).getName());
            }
        }
    }

    private void writeDataType(XMLWriter writer) {
        writer.writeOpenTag(MutationDeathTypeParser.MODEL_NAME,
                new Attribute.Default<String>(XMLParser.ID, DolloComponentOptions.DATA_NAME));
            writer.writeTag(MutationDeathTypeParser.EXTANT,
                    new Attribute.Default<String>(MutationDeathTypeParser.CODE, "1"), true);
            writer.writeTag(MutationDeathTypeParser.STATE,
                    new Attribute.Default<String>(MutationDeathTypeParser.CODE, "0"), true);
            writer.writeTag(MutationDeathTypeParser.AMBIGUITY,
                    new Attribute.Default<String>(MutationDeathTypeParser.CODE, "-"), true);  // Some users use this specification
            writer.writeTag(MutationDeathTypeParser.AMBIGUITY,
                    new Attribute.Default<String>(MutationDeathTypeParser.CODE, "?"), true);
        writer.writeCloseTag(MutationDeathTypeParser.MODEL_NAME);
    }


    private boolean includeStochasticDollo() {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model.isDolloModel()) {
                return true;
            }
        }
        return false;
    }

    @Override
    protected String getCommentLabel() {
        return "Stochastic Dollo";
    }

    private void writeDolloSubstitutionModels(XMLWriter writer, DolloComponentOptions component) {

        // generate tree likelihoods for stochastic Dollo partitions
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model.isDolloModel()) {
                writeDolloSubstitutionModel(partition, writer, component);
                writeDolloSiteModel(partition, writer, component);
                writeDolloTreeLikelihood(TreeLikelihoodParser.TREE_LIKELIHOOD, -1, partition, writer);
            }
        }
    }

    private void writeDolloSubstitutionModel(AbstractPartitionData partition, XMLWriter writer, DolloComponentOptions component) {
        String prefix = partition.getName() + ".";
//        String prefix = partition.getPrefix(); // TODO Fix
        writer.writeOpenTag(MutationDeathModelParser.MD_MODEL,
                new Attribute.Default<String>(XMLParser.ID, prefix + DolloComponentOptions.MODEL_NAME ));
        writeParameter(prefix + DolloComponentOptions.DEATH_RATE,
                component.getOptions().getParameter(prefix + DolloComponentOptions.DEATH_RATE), writer);
        writer.writeTag(DataType.DATA_TYPE, new Attribute.Default<String>(XMLParser.IDREF, DolloComponentOptions.DATA_NAME), true);
        writer.writeCloseTag(MutationDeathModelParser.MD_MODEL);
        writer.write("\n");
    }

    private void writeDolloSiteModel(AbstractPartitionData partition, XMLWriter writer, DolloComponentOptions components) {
        String prefix = partition.getName() + ".";
        writer.writeOpenTag(ALSSiteModelParser.ALS_SITE_MODEL,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, prefix + SiteModel.SITE_MODEL)});


        writer.writeOpenTag(GammaSiteModelParser.SUBSTITUTION_MODEL);
        writer.writeIDref(MutationDeathModelParser.MD_MODEL, prefix + DolloComponentOptions.MODEL_NAME);
        writer.writeCloseTag(GammaSiteModelParser.SUBSTITUTION_MODEL);

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
        if (model.hasCodonPartitions()) {
            writeParameter(GammaSiteModelParser.RELATIVE_RATE, "mu", model, writer);
        }

        if (model.isGammaHetero()) {
            writer.writeOpenTag(GammaSiteModelParser.GAMMA_SHAPE,
                    new Attribute.Default<String>(GammaSiteModelParser.GAMMA_CATEGORIES, "" + model.getGammaCategories()));
            writeParameter(prefix + "alpha", model, writer);
            writer.writeCloseTag(GammaSiteModelParser.GAMMA_SHAPE);
        }

        if (model.isInvarHetero()) {
            writeParameter(GammaSiteModelParser.PROPORTION_INVARIANT, "pInv", model, writer);
        }

        writer.writeCloseTag(ALSSiteModelParser.ALS_SITE_MODEL);
    }

//    private void writeDolloTreeLikelihoods(XMLWriter writer, DolloComponentOptions component) {
//
//        // generate tree likelihoods for stochastic Dollo partitions
//        for (AbstractPartitionData partition : options.dataPartitions) {
//            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
//            if (model.isDolloModel()) {
//                writeDolloTreeLikelihood(ALSTreeLikelihoodParser.LIKE_NAME, -1, partition, writer);
//            }
//        }
//    }

    private void writeDolloTreeLikelihood(String id, int num, AbstractPartitionData partition, XMLWriter writer) {
        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        String prefix = partition.getName() + ".";
        String oldPrefix = partition.getPrefix(); // TODO Get working

        writer.writeComment("Likelihood for tree given a stochastic Dollo model");

        writer.writeOpenTag(
                ALSTreeLikelihoodParser.LIKE_NAME,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, oldPrefix + id),
                        new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, true),
                        new Attribute.Default<Boolean>(ALSTreeLikelihoodParser.INTEGRATE_GAIN_RATE, true)}
        );

        if (!options.samplePriorOnly) {
            writer.writeIDref(SitePatternsParser.PATTERNS, partition.getPrefix() + SitePatternsParser.PATTERNS);
        } else {
            // We just need to use the dummy alignment
            if (partition instanceof PartitionData) {
                writer.writeIDref(AlignmentParser.ALIGNMENT, ((PartitionData) partition).getAlignment().getId());
            }
        }

        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
        writer.writeIDref(GammaSiteModel.SITE_MODEL, prefix + SiteModel.SITE_MODEL);

        writer.writeTag(ALSTreeLikelihoodParser.OBSERVATION_PROCESS,
                new Attribute.Default<String>(ALSTreeLikelihoodParser.OBSERVATION_TYPE,ALSTreeLikelihoodParser.ANY_TIP),
                true);

        switch (clockModel.getClockType()) {
            case STRICT_CLOCK:
                writer.writeIDref(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES, clockModel.getPrefix()
                        + BranchRateModel.BRANCH_RATES);
                break;
            case UNCORRELATED:
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES, clockModel.getPrefix()
                        + BranchRateModel.BRANCH_RATES);
                break;
            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(RandomLocalClockModelParser.LOCAL_BRANCH_RATES, clockModel.getPrefix()
                        + BranchRateModel.BRANCH_RATES);
                break;

            case AUTOCORRELATED:
                throw new UnsupportedOperationException("Autocorrelated relaxed clock model not implemented yet");
            default:
                throw new IllegalArgumentException("Unknown clock model");
        }
        writer.writeCloseTag(ALSTreeLikelihoodParser.LIKE_NAME);
    }

}// END: class