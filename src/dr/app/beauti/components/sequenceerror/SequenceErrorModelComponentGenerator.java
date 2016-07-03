/*
 * SequenceErrorModelComponentGenerator.java
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

package dr.app.beauti.components.sequenceerror;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.types.SequenceErrorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.alignment.HypermutantAlignment;
import dr.evomodel.tipstatesmodel.HypermutantErrorModel;
import dr.evomodelxml.tipstatesmodel.SequenceErrorModelParser;
import dr.evoxml.HypermutantAlignmentParser;
import dr.inference.model.ParameterParser;
import dr.inference.model.StatisticParser;
import dr.inferencexml.model.SumStatisticParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentGenerator extends BaseComponentGenerator {

    SequenceErrorModelComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        SequenceErrorModelComponentOptions comp = (SequenceErrorModelComponentOptions) options.getComponentOptions(SequenceErrorModelComponentOptions.class);


        if (!comp.usingSequenceErrorModel()) {
            return false;
        }

        switch (point) {
            case AFTER_PATTERNS:
            case AFTER_SITE_MODEL:
            case IN_TREE_LIKELIHOOD:
            case IN_FILE_LOG_PARAMETERS:
                return true;
        }
        return false;
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        SequenceErrorModelComponentOptions component = (SequenceErrorModelComponentOptions) options.getComponentOptions(SequenceErrorModelComponentOptions.class);

        switch (point) {
            case AFTER_PATTERNS:
                writeHypermutationAlignments(writer, component);
                break;
            case AFTER_SITE_MODEL:
                writeErrorModels(writer, component);
                break;
            case IN_TREE_LIKELIHOOD:
                AbstractPartitionData partition = (AbstractPartitionData) item;
                SequenceErrorType errorType = component.getSequenceErrorType(partition);
                if (errorType != SequenceErrorType.NO_ERROR) {
                    writer.writeIDref(SequenceErrorModelParser.SEQUENCE_ERROR_MODEL, partition.getPrefix() + "errorModel");
                }
                break;
            case IN_FILE_LOG_PARAMETERS:
                writeLogParameters(writer, component);
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    protected String getCommentLabel() {
        return "Sequence Error Model";
    }

    private void writeHypermutationAlignments(XMLWriter writer, SequenceErrorModelComponentOptions component) {
        for (AbstractPartitionData partition : options.getDataPartitions()) {
            String prefix = partition.getPrefix();//partition.getName() + ".";

            if (component.isHypermutation(partition)) {
                SequenceErrorType errorType = component.getSequenceErrorType(partition);

                final String errorTypeName;
                switch (errorType) {
                    case HYPERMUTATION_ALL:
                        errorTypeName = HypermutantAlignment.APOBECType.ALL.toString();
                        break;
                    case HYPERMUTATION_BOTH:
                        errorTypeName = HypermutantAlignment.APOBECType.BOTH.toString();
                        break;
                    case HYPERMUTATION_HA3F:
                        errorTypeName = HypermutantAlignment.APOBECType.HA3F.toString();
                        break;
                    case HYPERMUTATION_HA3G:
                        errorTypeName = HypermutantAlignment.APOBECType.HA3G.toString();
                        break;
                    default:
                        throw new RuntimeException("Unknown ErrorModelType: " + errorType.toString());
                }
                writer.writeOpenTag(
                        HypermutantAlignmentParser.HYPERMUTANT_ALIGNMENT,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + "hypermutants"),
                                new Attribute.Default<String>("type", errorTypeName)
                        }
                );

                writer.writeIDref("alignment", partition.getTaxonList().getId());

                writer.writeCloseTag(HypermutantAlignmentParser.HYPERMUTANT_ALIGNMENT);
            }
        }
    }

    private void writeErrorModels(XMLWriter writer, SequenceErrorModelComponentOptions component) {
        for (AbstractPartitionData partition : options.getDataPartitions()) {
            String prefix = partition.getPrefix();//partition.getName() + ".";

            SequenceErrorType errorType = component.getSequenceErrorType(partition);
            if (component.isHypermutation(partition)) {
                writer.writeOpenTag(
                        HypermutantErrorModel.HYPERMUTANT_ERROR_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + SequenceErrorModelComponentOptions.ERROR_MODEL)
                        }
                );

                writer.writeIDref(HypermutantAlignmentParser.HYPERMUTANT_ALIGNMENT, prefix + "hypermutants");

                writeParameter(HypermutantErrorModel.HYPERMUTATION_RATE, prefix + SequenceErrorModelComponentOptions.HYPERMUTION_RATE_PARAMETER, 1, writer);
                writeParameter(HypermutantErrorModel.HYPERMUTATION_INDICATORS, prefix + SequenceErrorModelComponentOptions.HYPERMUTANT_INDICATOR_PARAMETER, 1, writer);

                writer.writeCloseTag(HypermutantErrorModel.HYPERMUTANT_ERROR_MODEL);

                writer.writeOpenTag(SumStatisticParser.SUM_STATISTIC, new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, prefix + SequenceErrorModelComponentOptions.HYPERMUTANT_COUNT_STATISTIC),
                        new Attribute.Default<Boolean>(SumStatisticParser.ELEMENTWISE, true)});
                writer.writeIDref(ParameterParser.PARAMETER, prefix + SequenceErrorModelComponentOptions.HYPERMUTANT_INDICATOR_PARAMETER);
                writer.writeCloseTag(SumStatisticParser.SUM_STATISTIC);
            } else if (errorType != SequenceErrorType.NO_ERROR) {
                final String errorTypeName = (errorType == SequenceErrorType.AGE_TRANSITIONS ||
                        errorType == SequenceErrorType.BASE_TRANSITIONS ?
                        "transitions" : "all");

                writer.writeOpenTag(
                        SequenceErrorModelParser.SEQUENCE_ERROR_MODEL,
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, prefix + SequenceErrorModelComponentOptions.ERROR_MODEL),
                                new Attribute.Default<String>("type", errorTypeName)
                        }
                );

                if (component.hasAgeDependentRate(partition)) {
                    writeParameter(SequenceErrorModelComponentOptions.AGE_RATE, prefix + SequenceErrorModelComponentOptions.AGE_RATE_PARAMETER, 1, writer);
                }
                if (component.hasBaseRate(partition)) {
                    writeParameter(SequenceErrorModelComponentOptions.BASE_RATE, prefix + SequenceErrorModelComponentOptions.BASE_RATE_PARAMETER, 1, writer);
                }

                writer.writeCloseTag(SequenceErrorModelParser.SEQUENCE_ERROR_MODEL);
            }
        }
    }

    private void writeLogParameters(final XMLWriter writer, final SequenceErrorModelComponentOptions component) {
        for (AbstractPartitionData partition : options.getDataPartitions()) {
            String prefix = partition.getPrefix();//partition.getName() + ".";

            SequenceErrorType errorType = component.getSequenceErrorType(partition);
            if (errorType != SequenceErrorType.NO_ERROR) {
                if (component.isHypermutation(partition)) {
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + SequenceErrorModelComponentOptions.HYPERMUTION_RATE_PARAMETER);
                    writer.writeIDref(StatisticParser.STATISTIC, prefix + SequenceErrorModelComponentOptions.HYPERMUTANT_COUNT_STATISTIC);
                    writer.writeOpenTag(StatisticParser.STATISTIC,
                            new Attribute.Default<String>("name", "isHypermutated"));
                    writer.writeIDref(HypermutantErrorModel.HYPERMUTANT_ERROR_MODEL,
                            prefix + SequenceErrorModelComponentOptions.ERROR_MODEL);
                    writer.writeCloseTag(StatisticParser.STATISTIC);
                }

                if (component.hasAgeDependentRate(partition)) {
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + SequenceErrorModelComponentOptions.AGE_RATE_PARAMETER);
                }
                if (component.hasBaseRate(partition)) {
                    writer.writeIDref(ParameterParser.PARAMETER, prefix + SequenceErrorModelComponentOptions.BASE_RATE_PARAMETER);
                }
            }
        }
    }
}
