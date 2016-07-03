/*
 * AvgPosteriorIBDReporter.java
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

package dr.oldevomodel.ibd;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.oldevomodel.substmodel.AbstractSubstitutionModel;
import dr.oldevomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.NodePosteriorTreeLikelihood;
import dr.inference.model.*;
import dr.xml.*;

/**
 * Package: dr.evomodel.ibd
 * Description:  Computes for each tip the expected number of other tips IBD to it given the tip labels,
 * averaged over the full length of the alignment
 * <p/>
 * <p/>
 * Created by
 * avaleks (alexander.alekseyenko@gmail.com)
 * Date: 04-Aug-2008
 * Time: 13:46:33
 */
public class AvgPosteriorIBDReporter extends AbstractModel implements TreeTraitProvider {

    protected double[] ibdweights;
    protected double[][] ibdForward;
    protected double[][] ibdBackward;
    protected double[] diag;
    protected boolean weightsKnown;
    protected HKY substitutionModel;
    protected TreeModel treeModel;
    protected BranchRateModel branchRateModel;
    protected Parameter mutationParameter;
    protected NodePosteriorTreeLikelihood likelihoodReporter;
    protected double[] probabilities;

    AvgPosteriorIBDReporter(NodePosteriorTreeLikelihood likelihoodReporter, Parameter mutationParameter, TreeModel treeModel, BranchRateModel branchRateModel, AbstractSubstitutionModel substitutionModel) {
        super("AvgPosteriorIBDReporter");
        this.substitutionModel = (HKY) substitutionModel;
        addModel(this.substitutionModel);
        this.treeModel = treeModel;
        addModel(this.treeModel);
        this.branchRateModel = branchRateModel;
        addModel(this.branchRateModel);
        this.mutationParameter = mutationParameter;
        addVariable(this.mutationParameter);
        this.likelihoodReporter = likelihoodReporter;
        this.probabilities = new double[substitutionModel.getStateCount() * substitutionModel.getStateCount()];
    }

    public void forwardIBD() {
        int numNodes = treeModel.getNodeCount();
        int stateCount = substitutionModel.getStateCount();
        getDiagonalRates(diag);
        int patternCount = likelihoodReporter.getPatternCount();
        for (int nodeId = 0; nodeId < numNodes; ++nodeId) {
            NodeRef node = treeModel.getNode(nodeId);
            NodeRef parent = treeModel.getParent(node);
            likelihoodReporter.getNodeMatrix(nodeId, probabilities);
            double[] posteriors = likelihoodReporter.getPosteriors(nodeId);
            if (parent == null) { // handle the root

            } else if (treeModel.isExternal(node)) { // Handle the tip
                double branchTime = branchRateModel.getBranchRate(treeModel, node) * (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node));

                for (int state = 0; state < stateCount; ++state) {
                    double enorm = Math.exp(-diag[state] * branchTime) / probabilities[state + state * stateCount];
                    for (int pattern = 0; pattern < patternCount; ++pattern) {
                        ibdForward[nodeId][pattern * stateCount + state] = posteriors[pattern * stateCount + state] * enorm;
                    }
                }
            } else { // Handle internal node
                double branchTime = branchRateModel.getBranchRate(treeModel, node) * (treeModel.getNodeHeight(parent) - treeModel.getNodeHeight(node));

                int childCount = treeModel.getChildCount(node);
                for (int state = 0; state < stateCount; ++state) {
                    double enorm = Math.exp(-diag[state] * branchTime) / probabilities[state + state * stateCount];
                    for (int pattern = 0; pattern < patternCount; ++pattern) {
                        ibdForward[nodeId][pattern * stateCount + state] = 0;

                        for (int child = 0; child < childCount; ++child) {
                            int childNodeId = treeModel.getChild(node, child).getNumber();

                            ibdForward[nodeId][pattern * stateCount + state] += ibdForward[childNodeId][pattern * stateCount + state];
                        }
                        ibdForward[nodeId][pattern * stateCount + state] *= posteriors[pattern * stateCount + state] * enorm;
                    }
                }
            }
        }
    }

    public void backwardIBD(NodeRef node) {
        int stateCount = substitutionModel.getStateCount();
        int patternCount = likelihoodReporter.getPatternCount();
        if (node == null) {
            node = treeModel.getRoot();
            int nodeId = node.getNumber();
            for (int i = 0; i < patternCount * stateCount; ++i) {
                ibdBackward[nodeId][i] = 0;
            }
        }
        getDiagonalRates(diag);
        int childCount = treeModel.getChildCount(node);
        int nodeId = node.getNumber();

        double[] posteriors = likelihoodReporter.getPosteriors(nodeId);
        for (int child = 0; child < childCount; ++child) {
            NodeRef childNode = treeModel.getChild(node, child);
            int childNodeId = childNode.getNumber();
            likelihoodReporter.getNodeMatrix(childNodeId, probabilities);
            double branchTime = branchRateModel.getBranchRate(treeModel, childNode) * (treeModel.getNodeHeight(node) - treeModel.getNodeHeight(childNode));
            for (int pattern = 0; pattern < patternCount; ++pattern) {
                for (int state = 0; state < stateCount; ++state) {
                    ibdBackward[childNodeId][pattern * stateCount + state] = ibdBackward[nodeId][pattern * stateCount + state];
                    for (int sibling = 0; sibling < childCount; ++sibling) {
                        if (sibling != child) {
                            int siblingId = treeModel.getChild(node, sibling).getNumber();
                            ibdBackward[childNodeId][pattern * stateCount + state] += ibdForward[siblingId][pattern * stateCount + state];
                        }
                    }
                    ibdBackward[childNodeId][pattern * stateCount + state] *= posteriors[pattern * stateCount + state] * Math.exp(-diag[state] * branchTime) / probabilities[state + state * stateCount];
                }
            }
        }
        for (int child = 0; child < childCount; ++child) {
            NodeRef childNode = treeModel.getChild(node, child);
            backwardIBD(childNode);
        }

    }

    public void expectedIBD() {
        int stateCount = substitutionModel.getStateCount();
        int nodeCount = treeModel.getNodeCount();
        int patternCount = likelihoodReporter.getPatternCount();
        if (ibdweights == null) {
            ibdweights = new double[treeModel.getExternalNodeCount()];
            ibdForward = new double[nodeCount][stateCount * patternCount];
            ibdBackward = new double[nodeCount][stateCount * patternCount];
            diag = new double[stateCount];
        }


        forwardIBD();
        backwardIBD(null);
        int numTips = treeModel.getExternalNodeCount();
        double[] patternWeights = likelihoodReporter.getPatternWeights();
        double total = 0.0;
        for (int i = 0; i < patternCount; ++i) {
            total += patternWeights[i];
        }

        for (int i = 0; i < numTips; ++i) {
            double[] posteriors = likelihoodReporter.getPosteriors(i);

            ibdweights[i] = 0;
            for (int pattern = 0; pattern < patternCount; ++pattern) {
                for (int j = 0; j < stateCount; ++j) {
                    ibdweights[i] += ibdBackward[i][pattern * stateCount + j] * posteriors[pattern * stateCount + j] * patternWeights[pattern] / total;
                }
            }
        }
    }

    protected void getDiagonalRates(double[] diagonalRates) {
        double kappa = substitutionModel.getKappa();
        double[] freq = substitutionModel.getFrequencyModel().getFrequencies();
        double mutationRate = mutationParameter.getParameterValue(0);
        double beta = 0.5 / ((freq[0] + freq[2]) * (freq[1] + freq[3]) + kappa * (freq[0] * freq[2] + freq[1] * freq[3]));

        diagonalRates[0] = ((freq[1] + freq[3]) + freq[2] * kappa) * mutationRate * beta;
        diagonalRates[1] = ((freq[0] + freq[2]) + freq[3] * kappa) * mutationRate * beta;
        diagonalRates[2] = ((freq[1] + freq[3]) + freq[0] * kappa) * mutationRate * beta;
        diagonalRates[3] = ((freq[0] + freq[2]) + freq[1] * kappa) * mutationRate * beta;
    }

    TreeTrait avgPosteriorIBDWeight = new TreeTrait.D() {
        public String getTraitName() {
            return "AvgPosteriorIBDWeight";
        }

        public Intent getIntent() {
            return Intent.NODE;
        }

        public Double getTrait(Tree tree, NodeRef node) {
            if (!weightsKnown) {
                expectedIBD();
                weightsKnown = true;
            }
            if (tree.isExternal(node)) {
                int nodeNum = node.getNumber();
                return ibdweights[nodeNum] + 1;
            }
            return null;
        }
    };

    public TreeTrait[] getTreeTraits() {
        return new TreeTrait[] { avgPosteriorIBDWeight };
    }

    public TreeTrait getTreeTrait(String key) {
        // ignore the key - it must be the one they wanted, no?
        return avgPosteriorIBDWeight;
    }

    /**
     * The XML parser
     */

    public static final String IBD_REPORTER_LIKELIHOOD = "avgPosteriorIBDReporter";
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return IBD_REPORTER_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            Parameter mutationParameter = (Parameter) xo.getChild(Parameter.class);
            AbstractSubstitutionModel substitutionModel =
                    (AbstractSubstitutionModel) xo.getChild(AbstractSubstitutionModel.class);

            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
            if (branchRateModel == null) {
                branchRateModel = new DefaultBranchRateModel();
            }

            NodePosteriorTreeLikelihood likelihoodReporter = (NodePosteriorTreeLikelihood) xo.getChild(NodePosteriorTreeLikelihood.class);

            return new AvgPosteriorIBDReporter(likelihoodReporter, mutationParameter, treeModel, branchRateModel, substitutionModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents a reporter for average expected number of tips ibd conditional on observed patterns.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(AbstractSubstitutionModel.class),
                new ElementRule(Parameter.class),
                new ElementRule(NodePosteriorTreeLikelihood.class)
        };
    };

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == branchRateModel || model == treeModel || model == substitutionModel || model == likelihoodReporter) {
            weightsKnown = false;
        } else {
            System.err.println("Weird call back to IBDReporter from " + model.getModelName());
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == mutationParameter) {
            weightsKnown = false;
        } else {
            System.err.println("Weird call back to IBDReporter from " + variable.getVariableName());
        }
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }
}