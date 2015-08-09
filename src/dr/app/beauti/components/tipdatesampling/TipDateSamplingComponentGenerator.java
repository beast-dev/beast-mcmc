/*
 * TipDateSamplingComponentGenerator.java
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

package dr.app.beauti.components.tipdatesampling;

import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.types.TipDateSamplingType;
import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evoxml.TaxonParser;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TipDateSamplingComponentGenerator extends BaseComponentGenerator {

    public TipDateSamplingComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        TipDateSamplingComponentOptions comp = (TipDateSamplingComponentOptions)options.getComponentOptions(TipDateSamplingComponentOptions.class);

        if (comp.tipDateSamplingType == TipDateSamplingType.NO_SAMPLING) {
            return false;
        }

        switch (point) {
            case IN_TREE_MODEL:
            case IN_FILE_LOG_PARAMETERS:
                return true;
            case AFTER_TREE_MODEL:
                return options.getPartitionTreeModels().size() > 1 || comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT;
            case IN_MCMC_PRIOR:
                return comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT;
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        TipDateSamplingComponentOptions comp = (TipDateSamplingComponentOptions)options.getComponentOptions(TipDateSamplingComponentOptions.class);

        TaxonList taxa = comp.getTaxonSet();

        switch (point) {
            case IN_TREE_MODEL: {
                writeLeafHeightParameters(writer, (PartitionTreeModel)item, taxa);
            } break;
            case AFTER_TREE_MODEL:
                if (options.getPartitionTreeModels().size() > 1) {
                    // we have multiple treeModels with some or all the same taxa - create a JointParameter for each...

                    writeJointParameters(writer, taxa);
                }

                if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
                    writer.writeOpenTag("compoundParameter",
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, "treeModel.tipDates"),
                            }
                    );
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        writer.writeIDref(ParameterParser.PARAMETER, "age(" + taxon.getId() + ")");
                    }

                    writer.writeCloseTag("compoundParameter");
                }
                break;
            case IN_MCMC_PRIOR:
                if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY ||
                        comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_PRECISION) {
                    // nothing to do - individual parameter priors are written automatically
                } else if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {

                }
                break;
            case IN_FILE_LOG_PARAMETERS:
                if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY ||
                        comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_PRECISION) {
                    for (int i = 0; i < taxa.getTaxonCount(); i++) {
                        Taxon taxon = taxa.getTaxon(i);
                        writer.writeIDref(ParameterParser.PARAMETER, "age(" + taxon.getId() + ")");
                    }
                } else if (comp.tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
                    writer.writeIDref(ParameterParser.PARAMETER, "treeModel.tipDates");
                }
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }

    }

    private void writeJointParameters(XMLWriter writer, TaxonList taxa) {
        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);

            Set<PartitionTreeModel> treeModels = new HashSet<PartitionTreeModel>();
            for (PartitionTreeModel treeModel : options.getPartitionTreeModels()) {
                for (AbstractPartitionData data : options.getDataPartitions(treeModel)) {
                    if (data.getTaxonList().asList().contains(taxon)) {
                        treeModels.add(treeModel);
                    }
                }
            }

            // if we are sampling within precisions then only include this leaf if precision > 0
            if (treeModels.size() > 0) {
                writer.writeOpenTag("jointParameter",
                        new Attribute[]{
                                new Attribute.Default<String>(XMLParser.ID, "age(" + taxon.getId() + ")")
                        }
                );

                for (PartitionTreeModel treeModel : treeModels) {
                    writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, treeModel.getPrefix() + "age(" + taxon.getId() + ")"), true);
                }

                writer.writeCloseTag("jointParameter");
            }
        }
    }


    private void writeLeafHeightParameters(XMLWriter writer, PartitionTreeModel item, TaxonList taxa) {
        // only include this taxon as a leaf height if it found in this partition.
        PartitionTreeModel treeModel = (PartitionTreeModel)item;

        Set<Taxon> taxonSet = new HashSet<Taxon>();
        for (AbstractPartitionData data : options.getDataPartitions(treeModel)) {
            if (data.getTaxonList() != null) {
                for (Taxon taxon : data.getTaxonList()) {
                    taxonSet.add(taxon);
                }
            }
        }

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);

            if (taxonSet.contains(taxon)) {
                // if we are sampling within precisions then only include this leaf if precision > 0

                writer.writeOpenTag("leafHeight",
                        new Attribute[]{
                                new Attribute.Default<String>(TaxonParser.TAXON, taxon.getId()),
                        }
                );
                writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, treeModel.getPrefix() + "age(" + taxon.getId() + ")"), true);
                writer.writeCloseTag("leafHeight");
            }
        }
    }

    protected String getCommentLabel() {
        return "Tip date sampling";
    }

}