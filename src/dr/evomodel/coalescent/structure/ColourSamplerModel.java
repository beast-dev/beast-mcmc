/*
 * ColourSamplerModel.java
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

package dr.evomodel.coalescent.structure;

import dr.evolution.alignment.Alignment;
import dr.evolution.colouring.*;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @version $Id: ColourSamplerModel.java,v 1.14 2006/09/11 09:33:01 gerton Exp $
 */

public class ColourSamplerModel extends AbstractModel implements TreeColouringProvider, ModelListener, StatisticList {

    public static final String COLOUR_SAMPLER_MODEL = "colourSamplerModel";
    public static final String STRUCTURED_SAMPLER = "structuredSampler";
    public static final String NODE_BIAS = "nodeBias";
    public static final String BRANCH_BIAS = "branchBias";
    public static final String SECOND_ITERATION = "secondIteration";

    public ColourSamplerModel(TreeModel treeModel, ColourSampler colourSampler, MigrationModel migrationModel, MetaPopulationModel metaPopulationModel) {

        super(COLOUR_SAMPLER_MODEL);

        this.treeModel = treeModel;
        addModel(treeModel);

        this.migrationModel = migrationModel;
        addModel(migrationModel);

        this.metaPopulationModel = metaPopulationModel;
        // This is unusual because we don't want to addModel(metaPopulationModel) -
        // this is because the population sizes are just used as a bias in the proposal
        // distribution - a change in population size shouldn't force a recolouring by
        // this mechanism.

        addStatistic(migrationEventStatistic);
        //addStatistic(debugMigrationEventStatistic);

        addStatistic(rootColourStatistic);

        this.colourSampler = colourSampler;

    }

    public final TreeColouring getTreeColouring(Tree tree) {
        return getTreeColouring();
    }

    public final DefaultTreeColouring getTreeColouring() {
        if (treeColouring == null) {
            sample();
        }
        return treeColouring;
    }

    public final int[] getLeafColourCounts() {
        return colourSampler.getLeafColourCounts();
    }

    /**
     * Returns treeColouring and ensures that a proposal probability density has been assigned to it.
     *
     * @return treeColouring
     */
    public final DefaultTreeColouring getTreeColouringWithProbability() {

        DefaultTreeColouring tc = getTreeColouring();

        if (tc.hasProbability()) {
            return tc;
        }

        // We have a colouring, but it hasn't got a valid proposal probability.  This happens when parameters influencing the
        // proposal distribution have changed, or when local re-colouring moves have occurred. Re-compute the proposal probability

        // TODO HACK HACK HACK -- only using modern day population size
        // NOTE: If this is changed, StructuredCoalescentSampler should change too!  (GAL)
        //double[] N = metaPopulationModel.getPopulationSizes(0);

        double logP = colourSampler.getProposalProbability(tc, treeModel, migrationModel.getMigrationMatrix(), metaPopulationModel);
        tc.setLogProbabilityDensity(logP);
        return tc;
    }

    public final void resample() {
        treeColouring = null;
    }

    /**
     * Keeps the colouring, but reset its proposal probability (in response to a change in parameters)
     */
    public final void invalidateProposalProbability() {
        treeColouring = new DefaultTreeColouring(treeColouring);
    }

    /**
     * Colours the tree probabilistically with the given migration rates
     */
    private void sample() {

        // TODO HACK HACK HACK -- only using modern day population size
        // NOTE: If this is changed, StructuredCoalescentSampler should change too!  (GAL)
        //double[] N = metaPopulationModel.getPopulationSizes(0);

        treeColouring = colourSampler.sampleTreeColouring(treeModel, migrationModel.getMigrationMatrix(), metaPopulationModel);

        // StructuredCoalescent sc = new StructuredCoalescent();
        // StructuredIntervalList list = new ColouredTreeIntervals(treeModel, treeColouring);
        // double logL = sc.calculateLogLikelihood(treeColouring, list, migrationModel.getMigrationMatrix(), N);
        // double logP = treeColouring.getLogProbabilityDensity();

        //System.out.println("sampled tree colouring: " + logP + ", " + logL + " diff=" + (logL-logP));

        // treeColouring.checkColouring();

        listenerHelper.fireModelChanged(this);
    }

    /**
     * This function should be called to store the state of the
     * entire model. This makes the model state invalid until either
     * an acceptModelState or restoreModelState is called.
     */
    public void storeState() {
        storedTreeColouring = treeColouring;
    }

    /**
     * This function should be called to restore the state of the entire model.
     */
    public void restoreState() {
        treeColouring = storedTreeColouring;
    }

    /**
     * This function should be called to accept the state of the entire model
     */
    public void acceptState() {
        // Do nothing?
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        resample();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private final Statistic migrationEventStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "migrationEvents";
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int dim) {
            TreeColouring colouring = getTreeColouringWithProbability();
            return colouring.getColourChangeCount();
        }

    };

    private final Statistic rootColourStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "rootColour";
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int dim) {
            TreeColouring colouring = getTreeColouringWithProbability();
            return colouring.getNodeColour(colouring.getTree().getRoot());
        }

    };

    /*  for debugging
    private Statistic debugMigrationEventStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "migrationEvents";
        }

        public int getDimension() { return 4; }

        public double getStatisticValue(int dim) {
            TreeColouring colouring = getTreeColouringWithProbability();

            if (dim == 0) return colouring.getColourChangeCount();

            if (dim == 1) return colouring.getLogProbabilityDensity();

            StructuredCoalescent sc = new StructuredCoalescent();

            StructuredIntervalList list = new ColouredTreeIntervals(treeModel, colouring);
            
            //double[] N = metaPopulationModel.getPopulationSizes(0);

            double logL = sc.calculateLogLikelihood(colouring, list, migrationModel.getMigrationMatrix(), metaPopulationModel);

            if (dim == 2) return logL;

            return logL-colouring.getLogProbabilityDensity();
        }

    };
 */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COLOUR_SAMPLER_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            MigrationModel migrationModel = (MigrationModel) xo.getChild(MigrationModel.class);

            MetaPopulationModel metaPopulationModel = (MetaPopulationModel) xo.getChild(MetaPopulationModel.class);

            boolean structuredSampler = xo.getAttribute(STRUCTURED_SAMPLER, true);
            boolean branchBias = xo.getAttribute(BRANCH_BIAS, true);
            boolean nodeBias = xo.getAttribute(NODE_BIAS, true);
            boolean secondIteration = xo.getAttribute(SECOND_ITERATION, false);

            ColourSampler colourSampler;
            if (xo.hasChildNamed("colours")) {
                XMLObject cxo = xo.getChild("colours");

                Taxa[] colourTaxa = new Taxa[cxo.getChildCount()];
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    colourTaxa[i] = (Taxa) cxo.getChild(i);
                }

                if (structuredSampler) {
                    colourSampler = new StructuredColourSampler(colourTaxa, treeModel, nodeBias, branchBias, secondIteration);
                } else {
                    colourSampler = new BasicColourSampler(colourTaxa, treeModel);
                }
            } else {
                Alignment alignment = (Alignment) xo.getChild(Alignment.class);
                if (structuredSampler) {
                    colourSampler = new StructuredColourSampler(alignment, treeModel, nodeBias, branchBias, secondIteration);
                } else {
                    colourSampler = new BasicColourSampler(alignment, treeModel);
                }
            }

            ColourSamplerModel colourSamplerModel = new ColourSamplerModel(treeModel, colourSampler, migrationModel, metaPopulationModel);

            if (structuredSampler) {
                Logger.getLogger("dr.evomodel").info("Creating colour sampler model with 2 colours");
                if (!nodeBias) {
                    Logger.getLogger("dr.evomodel").info(" Colour sampler has node biases switched off");
                }
                if (!branchBias) {
                    Logger.getLogger("dr.evomodel").info(" Colour sampler has branch biases switched off");
                }
            } else {
                Logger.getLogger("dr.evomodel").info("Creating basic 2-colour sampler");
            }

            return colourSamplerModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a likelihood function for transmission.";
        }

        public Class getReturnType() {
            return StructuredCoalescentLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new XORRule(
                        new ElementRule("colours", new XMLSyntaxRule[]{
                                new ElementRule(Taxa.class, "Taxa for each subsequent colour (after 0).", 1, Integer.MAX_VALUE),
                        }),
                        new ElementRule(Alignment.class, "The alignment.")),
                new ElementRule(TreeModel.class, "The tree."),
                new ElementRule(MigrationModel.class, "The migration model."),
                new ElementRule(MetaPopulationModel.class, "The metapopulation model,")
        };
    };


    private final MetaPopulationModel metaPopulationModel;
    private final MigrationModel migrationModel;
    private final ColourSampler colourSampler;
    private final TreeModel treeModel;

    private DefaultTreeColouring treeColouring = null;
    private DefaultTreeColouring storedTreeColouring = null;

}
