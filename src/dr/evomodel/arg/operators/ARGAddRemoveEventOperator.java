/*
 * ARGAddRemoveEventOperator.java
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

/*
 * AddRemoveSubtreeOperator.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.arg.operators;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGModel.Node;
import dr.evomodel.arg.ARGPartitionLikelihood;
import dr.evomodel.arg.ARGRatePrior;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.logging.Logger;


/**
 * Adds and removes re-assortment events.
 *
 * @author Erik Bloomquist
 * @author Marc Suchard
 */

public class ARGAddRemoveEventOperator extends AbstractCoercableOperator {

    public static final String ADD_PROBABILITY = "addProbability";
    public static final String ARG_EVENT_OPERATOR = "ARGEventOperator";
    public static final String INTERNAL_NODES = "internalNodes";
    public static final String INTERNAL_AND_ROOT = "internalNodesPlusRoot";
    public static final String NODE_RATES = TreeModelParser.NODE_RATES;
    public static final String BELOW_ROOT_PROBABILITY = "belowRootProbability";
    public static final String FLIP_MEAN = "flipMean";
    public static final String JOINT_PARTITIONING = "jointPartitioning";
    public static final String RELAXED = "relaxed";

    public static final double LOG_TWO = Math.log(2.0);


    private boolean relaxed = false;

    private ARGModel arg = null;
    private double size = 0.0;  //Translates into add probability of 50%
    private double probBelowRoot = 0.9; //Transformed in constructor for computational efficiency

    private ARGPartitionLikelihood partLike;
    private ARGRatePrior ratePrior;


    //	private int mode = CoercableMCMCOperator.COERCION_OFF;
    private CompoundParameter internalNodeParameters;
    private CompoundParameter internalAndRootNodeParameters;
    private CompoundParameter nodeRates;

    private int tossSize = -1;


    public ARGAddRemoveEventOperator(ARGModel arg, int weight, double size, CoercionMode mode,
                                     CompoundParameter param1,
                                     CompoundParameter param2,
                                     CompoundParameter param3,
                                     double belowRootProbability,
                                     ARGPartitionLikelihood partLike,
                                     ARGRatePrior ratePrior,
                                     int tossSize) {
        super(mode);
        this.arg = arg;
        this.size = size;
        this.internalNodeParameters = param1;
        this.internalAndRootNodeParameters = param2;
        this.nodeRates = param3;

        this.partLike = partLike;
        this.ratePrior = ratePrior;

        if (ratePrior != null) {
            relaxed = true;
        }

//        this.isRecombination = arg.isRecombinationPartitionType();
////		this.mode = mode;
//        this.flipMean = flipMean;
//
//        if (flipMean > -1) {
//            this.pd = new PoissonDistribution(flipMean);
//            logPoissonNormalizingFactor = Math.log(pd.cdf(arg.getNumberOfPartitions() - 1) - pd.cdf(0));
//        }
//
//        if (isRecombination) {
//            logPartitionSize = Math.log(arg.getNumberOfPartitions() - 1);
//        } else if (arg.getNumberOfPartitions() > 40) {
//            //The -1 factor only matters for smaller problems,
//            //so we ignore when partition size gets big.
//            //Difference is very, very tiny.
//            logPartitionSize = (arg.getNumberOfPartitions() - 1) * LOG_TWO;
//        } else {
//            logPartitionSize = Math.log(Math.pow(2, arg.getNumberOfPartitions() - 1) - 1);
//        }

        setWeight(weight);

        this.probBelowRoot = belowRootProbability;

        //This is for computational efficiency
        probBelowRoot = -Math.log(1 - Math.sqrt(probBelowRoot));

        this.tossSize = tossSize;
    }


    /**
     * Do a add/remove re-assortment node operation
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        double logq = 0;

        try {
            if (MathUtils.nextDouble() < 1.0 / (1 + Math.exp(-size))) {
                logq = AddOperation() - size;
            } else {
                logq = RemoveOperation() + size;
            }
        } catch (NoReassortmentEventException nree) {
            throw new RuntimeException("");
        }

        assert !Double.isInfinite(logq) && !Double.isNaN(logq);


        return logq;
    }

    private double AddOperation() {
        double logHastings = 0;


        double treeHeight = arg.getNodeHeight(arg.getRoot());
        double newBifurcationHeight = Double.POSITIVE_INFINITY;
        double newReassortmentHeight = Double.POSITIVE_INFINITY;

        double theta = probBelowRoot / treeHeight;

        while (newBifurcationHeight > treeHeight && newReassortmentHeight > treeHeight) {
            newBifurcationHeight = MathUtils.nextExponential(theta);
            newReassortmentHeight = MathUtils.nextExponential(theta);
        }

        logHastings += theta * (newBifurcationHeight + newReassortmentHeight) - LOG_TWO
                - 2.0 * Math.log(theta) + Math.log(1 - Math.exp(-2.0 * treeHeight * theta));

//		This is the ugly mixture proposal
//		double meanRoot = 4.0 / treeHeight;
//		double case1 = 0.95;
//
//		if(MathUtils.nextDouble() < case1){
//			newReassortmentHeight = MathUtils.nextDouble() * treeHeight;
//			newBifurcationHeight = MathUtils.nextDouble() * treeHeight;
//
//			logHastings += 2.0*Math.log(treeHeight) - Math.log(2.0*case1);
//		}else{
//			newReassortmentHeight = MathUtils.nextDouble() * treeHeight;
//			double additional = MathUtils.nextExponential(meanRoot);
//			logHastings += Math.log(treeHeight) + additional*meanRoot -
//				Math.log((1-case1)*meanRoot);
//
//			newBifurcationHeight = additional + treeHeight;
//		}

        if (newBifurcationHeight < newReassortmentHeight) {
            double temp = newBifurcationHeight;
            newBifurcationHeight = newReassortmentHeight;
            newReassortmentHeight = temp;
        }

        //2. Find the possible re-assortment and bifurcation points.
        ArrayList<NodeRef> potentialBifurcationChildren = new ArrayList<NodeRef>();
        ArrayList<NodeRef> potentialReassortmentChildren = new ArrayList<NodeRef>();

        int totalPotentialBifurcationChildren = findPotentialAttachmentPoints(
                newBifurcationHeight, potentialBifurcationChildren);
        int totalPotentialReassortmentChildren = findPotentialAttachmentPoints(
                newReassortmentHeight, potentialReassortmentChildren);

        assert totalPotentialBifurcationChildren > 0;
        assert totalPotentialReassortmentChildren > 0;

        logHastings += Math.log((double) potentialBifurcationChildren.size() *
                potentialReassortmentChildren.size());

        //3.  Choose your re-assortment location.
        Node reassortChild = (Node) potentialReassortmentChildren.get(MathUtils
                .nextInt(totalPotentialReassortmentChildren));

        Node reassortLeftParent = reassortChild.leftParent;
        Node reassortRightParent = reassortChild.rightParent;
        Node reassortSplitParent = reassortChild.leftParent;

        boolean splitReassortLeftParent = true;

        if (!reassortChild.bifurcation) {
            boolean[] tester = {arg.getNodeHeight(reassortLeftParent) > newReassortmentHeight,
                    arg.getNodeHeight(reassortRightParent) > newReassortmentHeight};

            if (tester[0] && tester[1]) {
                if (MathUtils.nextBoolean()) {
                    splitReassortLeftParent = false;
                    reassortSplitParent = reassortRightParent;
                }
//        		 logHastings += LOG_TWO;
            } else if (tester[0]) {
                //nothing to do, setup above
            } else if (tester[1]) {
                splitReassortLeftParent = false;
                reassortSplitParent = reassortRightParent;
            } else {
                assert false;
            }
        }


//        if (recParentL != recParentR) {
//            boolean[] tester = {arg.getNodeHeight(recParentL) > newReassortmentHeight,
//                    arg.getNodeHeight(recParentR) > newReassortmentHeight};
//
//            if (tester[0] && tester[1]) {
//                if (MathUtils.nextBoolean()) {
//                    recParent = recParentR;
//                }
//
//                logHastings += LOG_TWO;
//            } else if (tester[0]) {
//                recParent = recParentL;
//            } else {
//                recParent = recParentR;
//            }
//        }

        //4. Choose your bifurcation location.

        Node bifurcationChild = (Node) potentialBifurcationChildren.get(MathUtils
                .nextInt(potentialBifurcationChildren.size()));


        Node bifurcationLeftParent = bifurcationChild.leftParent;
        Node bifurcationRightParent = bifurcationChild.rightParent;
        Node bifurcationSplitParent = bifurcationLeftParent;

        boolean splitBifurcationLeftParent = true;

        if (!bifurcationChild.bifurcation) {
            boolean[] tester = {arg.getNodeHeight(bifurcationLeftParent) > newBifurcationHeight,
                    arg.getNodeHeight(bifurcationRightParent) > newBifurcationHeight};

            if (tester[0] && tester[1]) {
                if (MathUtils.nextBoolean()) {
                    splitBifurcationLeftParent = false;
                    bifurcationSplitParent = bifurcationRightParent;
                }
//    		  logHastings += LOG_TWO;
            } else if (tester[0]) {
                //nothing to do
            } else if (tester[1]) {
                splitBifurcationLeftParent = false;
                bifurcationSplitParent = bifurcationRightParent;
            } else {
                assert false;
            }
        }

        boolean attachNewReassortNewBifurcationThroughLeft = MathUtils.nextBoolean();

        //During the delete step, this guy gets cancelled out!
        logHastings += LOG_TWO;


//        if (sisParentL != sisParentR) {
//            boolean[] tester = {arg.getNodeHeight(sisParentL) > newBifurcationHeight,
//                    arg.getNodeHeight(sisParentR) > newBifurcationHeight};
//
//            if (tester[0] && tester[1]) {
//                if (MathUtils.nextBoolean()) {
//                    sisParent = sisParentR;
//                }
//                logHastings += LOG_TWO;
//            } else if (tester[0]) {
//                sisParent = sisParentL;
//            } else {
//                sisParent = sisParentR;
//            }
//        }

        //5. Make the new nodes.
        //Note: The height stuff is taken care of below.

        Node newReassortment = arg.new Node();
        newReassortment.bifurcation = false;

        double[] reassortmentValues = {1.0};

        if (relaxed) {
            reassortmentValues = ratePrior.generateValues();
        }

//        logHastings += ratePrior.getAddHastingsRatio(reassortmentValues);

        newReassortment.rateParameter = new Parameter.Default(reassortmentValues);

        newReassortment.rateParameter.addBounds(new Parameter.DefaultBounds(
                Double.POSITIVE_INFINITY, 0, reassortmentValues.length));

        newReassortment.number = arg.getNodeCount() + 1;


        Node newBifurcation = arg.new Node();

        double[] bifurcationValues = {1.0};

        if (relaxed) {
            bifurcationValues = ratePrior.generateValues();
            logHastings += ratePrior.getAddHastingsRatio(bifurcationValues);
        }

        newBifurcation.rateParameter = new Parameter.Default(bifurcationValues);

        newBifurcation.rateParameter.addBounds(new Parameter.DefaultBounds(
                Double.POSITIVE_INFINITY, 0, bifurcationValues.length));

        newBifurcation.number = arg.getNodeCount();

        //6. Begin editing the tree.
        arg.beginTreeEdit();

        adjustRandomPartitioning();

        if (newBifurcationHeight < treeHeight) {

            if (bifurcationChild == reassortChild) {
                if (bifurcationChild.bifurcation) {
                    bifurcationChild.leftParent = bifurcationChild.rightParent = newReassortment;
                    newReassortment.leftChild = newReassortment.rightChild = bifurcationChild;
                    newReassortment.leftParent = newReassortment.rightParent = newBifurcation;
                    newBifurcation.leftChild = newBifurcation.rightChild = newReassortment;
                    newBifurcation.leftParent = newBifurcation.rightParent = bifurcationSplitParent;

                    if (bifurcationSplitParent.bifurcation) {
                        if (bifurcationSplitParent.leftChild == bifurcationChild) {
                            bifurcationSplitParent.leftChild = newBifurcation;
                        } else {
                            bifurcationSplitParent.rightChild = newBifurcation;
                        }
                    } else {
                        bifurcationSplitParent.leftChild = bifurcationSplitParent.rightChild = newBifurcation;
                    }

                    logHastings -= LOG_TWO;
                } else {
                    if (splitBifurcationLeftParent && splitReassortLeftParent) {
                        bifurcationChild.leftParent = newReassortment;
                        newReassortment.leftChild = newReassortment.rightChild = bifurcationChild;
                        newReassortment.leftParent = newReassortment.rightParent = newBifurcation;

                        newBifurcation.leftChild = newBifurcation.rightChild = newReassortment;
                        newBifurcation.leftParent = newBifurcation.rightParent = bifurcationSplitParent;

                        if (bifurcationSplitParent.bifurcation) {
                            if (bifurcationSplitParent.leftChild == bifurcationChild) {
                                bifurcationSplitParent.leftChild = newBifurcation;
                            } else {
                                bifurcationSplitParent.rightChild = newBifurcation;
                            }
                        } else {
                            bifurcationSplitParent.leftChild = bifurcationSplitParent.rightChild = newBifurcation;
                        }


                    } else if (splitBifurcationLeftParent) {
                        //bifurcation on left, reassortment on right

                        bifurcationChild.leftParent = newBifurcation;
                        bifurcationChild.rightParent = newReassortment;

                        newBifurcation.leftChild = bifurcationChild;
                        newBifurcation.rightChild = newReassortment;
                        newBifurcation.leftParent = newBifurcation.rightParent = bifurcationSplitParent;

                        if (bifurcationSplitParent.bifurcation) {
                            if (bifurcationSplitParent.leftChild == bifurcationChild) {
                                bifurcationSplitParent.leftChild = newBifurcation;
                            } else {
                                bifurcationSplitParent.rightChild = newBifurcation;
                            }
                        } else {
                            bifurcationSplitParent.leftChild = bifurcationSplitParent.rightChild = newBifurcation;
                        }

                        newReassortment.leftChild = newReassortment.rightChild = bifurcationChild;

                        if (attachNewReassortNewBifurcationThroughLeft) {
                            newReassortment.leftParent = newBifurcation;
                            newReassortment.rightParent = reassortSplitParent;
                        } else {
                            newReassortment.rightParent = newBifurcation;
                            newReassortment.leftParent = reassortSplitParent;
                        }

                        if (reassortSplitParent.bifurcation) {
                            if (reassortSplitParent.leftChild == reassortChild) {
                                reassortSplitParent.leftChild = newReassortment;
                            } else {
                                reassortSplitParent.rightChild = newReassortment;
                            }
                        } else {
                            reassortSplitParent.leftChild = reassortSplitParent.rightChild = newReassortment;
                        }


                    } else if (splitReassortLeftParent) {
                        //bifurcation on right, reassortment on left

                        bifurcationChild.rightParent = newBifurcation;
                        bifurcationChild.leftParent = newReassortment;

                        newBifurcation.leftChild = bifurcationChild;
                        newBifurcation.rightChild = newReassortment;
                        newBifurcation.leftParent = newBifurcation.rightParent = bifurcationSplitParent;

                        if (bifurcationSplitParent.bifurcation) {
                            if (bifurcationSplitParent.leftChild == bifurcationChild) {
                                bifurcationSplitParent.leftChild = newBifurcation;
                            } else {
                                bifurcationSplitParent.rightChild = newBifurcation;
                            }
                        } else {
                            bifurcationSplitParent.leftChild = bifurcationSplitParent.rightChild = newBifurcation;
                        }

                        newReassortment.leftChild = newReassortment.rightChild = bifurcationChild;

                        if (attachNewReassortNewBifurcationThroughLeft) {
                            newReassortment.leftParent = newBifurcation;
                            newReassortment.rightParent = reassortSplitParent;
                        } else {
                            newReassortment.rightParent = newBifurcation;
                            newReassortment.leftParent = reassortSplitParent;
                        }

                        if (reassortSplitParent.bifurcation) {
                            if (reassortSplitParent.leftChild == reassortChild) {
                                reassortSplitParent.leftChild = newReassortment;
                            } else {
                                reassortSplitParent.rightChild = newReassortment;
                            }
                        } else {
                            reassortSplitParent.leftChild = reassortSplitParent.rightChild = newReassortment;
                        }


                    } else {

                        bifurcationChild.rightParent = newReassortment;
                        newReassortment.leftChild = newReassortment.rightChild = bifurcationChild;
                        newReassortment.leftParent = newReassortment.rightParent = newBifurcation;

                        newBifurcation.leftChild = newBifurcation.rightChild = newReassortment;
                        newBifurcation.leftParent = newBifurcation.rightParent = bifurcationSplitParent;

                        if (bifurcationSplitParent.bifurcation) {
                            if (bifurcationSplitParent.leftChild == bifurcationChild) {
                                bifurcationSplitParent.leftChild = newBifurcation;
                            } else {
                                bifurcationSplitParent.rightChild = newBifurcation;
                            }
                        } else {
                            bifurcationSplitParent.leftChild = bifurcationSplitParent.rightChild = newBifurcation;
                        }

                        logHastings -= LOG_TWO;

                    }
                }

            } else {


                newReassortment.leftChild = newReassortment.rightChild = reassortChild;
                newBifurcation.leftParent = newBifurcation.rightParent = bifurcationSplitParent;
                newBifurcation.leftChild = newReassortment;
                newBifurcation.rightChild = bifurcationChild;

                if (attachNewReassortNewBifurcationThroughLeft) {
                    newReassortment.leftParent = newBifurcation;
                    newReassortment.rightParent = reassortSplitParent;
                } else {
                    newReassortment.rightParent = newBifurcation;
                    newReassortment.leftParent = reassortSplitParent;
                }

                if (reassortChild.bifurcation) {
                    reassortChild.leftParent = reassortChild.rightParent = newReassortment;
                } else if (splitReassortLeftParent) {
                    reassortChild.leftParent = newReassortment;
                } else {
                    reassortChild.rightParent = newReassortment;
                }

                if (reassortSplitParent.bifurcation) {
                    if (reassortSplitParent.leftChild == reassortChild) {
                        reassortSplitParent.leftChild = newReassortment;
                    } else {
                        reassortSplitParent.rightChild = newReassortment;
                    }
                } else {
                    reassortSplitParent.leftChild = reassortSplitParent.rightChild = newReassortment;
                }

                if (bifurcationChild.bifurcation) {
                    bifurcationChild.leftParent = bifurcationChild.rightParent = newBifurcation;
                } else if (splitBifurcationLeftParent) {
                    bifurcationChild.leftParent = newBifurcation;
                } else {
                    bifurcationChild.rightParent = newBifurcation;
                }

                if (bifurcationSplitParent.bifurcation) {
                    if (bifurcationSplitParent.leftChild == bifurcationChild) {
                        bifurcationSplitParent.leftChild = newBifurcation;
                    } else {
                        bifurcationSplitParent.rightChild = newBifurcation;
                    }
                } else {
                    bifurcationSplitParent.leftChild = bifurcationSplitParent.rightChild = newBifurcation;
                }


            }

            Parameter partition = new Parameter.Default(arg.getNumberOfPartitions());
            drawRandomPartitioning(partition);


            newReassortment.partitioning = partition;

            newBifurcation.heightParameter = new Parameter.Default(newBifurcationHeight);
            newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
            newBifurcation.setupHeightBounds();
            newReassortment.setupHeightBounds();


            arg.expandARG(newBifurcation, newReassortment,
                    internalNodeParameters,
                    internalAndRootNodeParameters,
                    nodeRates);


//             arg.expandARGWithRecombinant(newBifurcation, newReassortment,
//                     internalNodeParameters,
//                     internalAndRootNodeParameters,
//                     nodeRates);
            assert nodeCheck() : arg.toARGSummary();

        } else {

            assert newReassortmentHeight < treeHeight;


            //New bifurcation takes the place of the old root.
            //Much easier to program.

            newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
            newReassortment.setupHeightBounds();

            bifurcationChild = newBifurcation;
            if (arg.isRoot(reassortSplitParent))
                reassortSplitParent = newBifurcation;


            Node root = (Node) arg.getRoot();
            Node rootLeftChild = root.leftChild;
            Node rootRightChild = root.rightChild;

            arg.singleRemoveChild(root, rootLeftChild);
            arg.singleRemoveChild(root, rootRightChild);
            arg.singleAddChild(newBifurcation, rootLeftChild);
            arg.singleAddChild(newBifurcation, rootRightChild);

            if (reassortSplitParent.isBifurcation())
                arg.singleRemoveChild(reassortSplitParent, reassortChild);
            else
                arg.doubleRemoveChild(reassortSplitParent, reassortChild);

            arg.doubleAddChild(newReassortment, reassortChild);
            arg.singleAddChild(root, newBifurcation);

            Parameter partitioning = new Parameter.Default(arg.getNumberOfPartitions());
            drawRandomPartitioning(partitioning);


            arg.addChildAsRecombinant(root, reassortSplitParent, newReassortment, partitioning);

            if (attachNewReassortNewBifurcationThroughLeft) {
                newReassortment.leftParent = root;
                newReassortment.rightParent = reassortSplitParent;
            } else {
                newReassortment.leftParent = reassortSplitParent;
                newReassortment.rightParent = root;
            }

            newBifurcation.heightParameter = new Parameter.Default(root.getHeight());

            newBifurcation.setupHeightBounds();
            root.heightParameter.setParameterValue(0, newBifurcationHeight);


            arg.expandARG(newBifurcation, newReassortment,
                    internalNodeParameters,
                    internalAndRootNodeParameters,
                    nodeRates);


//         arg.expandARGWithRecombinant(newBifurcation, newReassortment,
//                 internalNodeParameters,
//                 internalAndRootNodeParameters,
//                 nodeRates);

            assert nodeCheck();


        }


        //6a. This is when we do not create a new root.
//        if (newBifurcationHeight < treeHeight) {
//            newBifurcation.heightParameter = new Parameter.Default(newBifurcationHeight);
//            newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
//            newBifurcation.setupHeightBounds();
//            newReassortment.setupHeightBounds();
//
//            if (sisParent.bifurcation)
//                arg.singleRemoveChild(sisParent, bifurcationChild);
//            else
//                arg.doubleRemoveChild(sisParent, bifurcationChild);
//            if (bifurcationChild != reassortChild) {
//                if (recParent.bifurcation)
//                    arg.singleRemoveChild(recParent, reassortChild);
//                else
//                    arg.doubleRemoveChild(recParent, reassortChild);
//            }
//            if (sisParent.bifurcation)
//                arg.singleAddChild(sisParent, newBifurcation);
//            else
//                arg.doubleAddChild(sisParent, newBifurcation);
//            if (bifurcationChild != reassortChild)
//                arg.singleAddChild(newBifurcation, bifurcationChild);
//            arg.doubleAddChild(newReassortment, reassortChild);
//
//            Parameter partitioning = new Parameter.Default(arg.getNumberOfPartitions());
//            drawRandomPartitioning(partitioning);
//
//            if (bifurcationChild != reassortChild) {
//                arg.addChildAsRecombinant(newBifurcation, recParent,
//                        newReassortment, partitioning);
//            } else {
//                arg.addChildAsRecombinant(newBifurcation, newBifurcation,
//                        newReassortment, partitioning);
//            }
//            arg.expandARGWithRecombinant(newBifurcation, newReassortment,
//                    internalNodeParameters,
//                    internalAndRootNodeParameters,
//                    nodeRates);
//            assert nodeCheck();
//
//            //6b. But here we do.
//        } else if (newReassortmentHeight < treeHeight) {
//
//            newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
//            newReassortment.setupHeightBounds();
//
//            bifurcationChild = newBifurcation;
//            if (arg.isRoot(recParent))
//                recParent = newBifurcation;
//
//
//            Node root = (Node) arg.getRoot();
//            Node rootLeftChild = root.leftChild;
//            Node rootRightChild = root.rightChild;
//
//            arg.singleRemoveChild(root, rootLeftChild);
//            arg.singleRemoveChild(root, rootRightChild);
//            arg.singleAddChild(newBifurcation, rootLeftChild);
//            arg.singleAddChild(newBifurcation, rootRightChild);
//
//            if (recParent.isBifurcation())
//                arg.singleRemoveChild(recParent, reassortChild);
//            else
//                arg.doubleRemoveChild(recParent, reassortChild);
//
//            arg.doubleAddChild(newReassortment, reassortChild);
//            arg.singleAddChild(root, newBifurcation);
//
//            Parameter partitioning = new Parameter.Default(arg.getNumberOfPartitions());
//            drawRandomPartitioning(partitioning);
//
//
//            arg.addChildAsRecombinant(root, recParent, newReassortment, partitioning);
//
//            newBifurcation.heightParameter = new Parameter.Default(root.getHeight());
//
//            newBifurcation.setupHeightBounds();
//            root.heightParameter.setParameterValue(0, newBifurcationHeight);
//
//
//            arg.expandARGWithRecombinant(newBifurcation, newReassortment,
//                    internalNodeParameters, internalAndRootNodeParameters,
//                    nodeRates);
//
//            assert nodeCheck();
//
//        } else {
//
//            Node root = (Node) arg.getRoot();
//            Node rootLeftChild = root.leftChild;
//            Node rootRightChild = root.rightChild;
//
//            arg.singleRemoveChild(root, rootLeftChild);
//            arg.singleRemoveChild(root, rootRightChild);
//            arg.singleAddChild(newBifurcation, rootLeftChild);
//            arg.singleAddChild(newBifurcation, rootRightChild);
//
//            arg.doubleAddChild(newReassortment, newBifurcation);
//            arg.doubleAddChild(root, newReassortment);
//
//            Parameter partitioning = new Parameter.Default(arg.getNumberOfPartitions());
//            drawRandomPartitioning(partitioning);
//
//            newReassortment.partitioning = partitioning;
//
//            newBifurcation.heightParameter = new Parameter.Default(arg.getNodeHeight(root));
//            newReassortment.heightParameter = new Parameter.Default(newReassortmentHeight);
//            root.heightParameter.setParameterValueQuietly(0, newBifurcationHeight);
//
//            newBifurcation.setupHeightBounds();
//            newReassortment.setupHeightBounds();
//
//            arg.expandARGWithRecombinant(newBifurcation, newReassortment,
//                    internalNodeParameters, internalAndRootNodeParameters,
//                    nodeRates);
//
//            assert nodeCheck();
//
//        }

        arg.pushTreeSizeIncreasedEvent();

        arg.endTreeEdit();
        try {
            arg.checkTreeIsValid();
        } catch (MutableTree.InvalidTreeException ite) {
            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                    + "\n" + TreeUtils.uniqueNewick(arg, arg.getRoot()));
        }

        //Do all the backwards stuff now. :(

        assert nodeCheck();

        logHastings -= Math.log((double) findPotentialNodesToRemove(null));
//		logHastings -= Math.log((double)arg.getReassortmentNodeCount());

//		if (newReassortment.leftParent != newReassortment.rightParent){
//			if(newReassortment.leftParent.bifurcation
//				&& newReassortment.rightParent.bifurcation)
//				logHastings -= LOG_TWO;
//		}

        assert nodeCheck();

        assert !Double.isNaN(logHastings) && !Double.isInfinite(logHastings);

        if (newReassortment.leftParent.bifurcation && newReassortment.rightParent.bifurcation
                && newReassortment.leftParent != newReassortment.rightParent) {
            logHastings -= LOG_TWO;
        }

        //You're done, return the hastings ratio!

//		System.out.println(logHastings);

        logHastings += getPartitionAddHastingsRatio(newReassortment.partitioning.getParameterValues());

        return logHastings;
    }

    private int findPotentialAttachmentPoints(double time, ArrayList<NodeRef> list) {
        int count = 0;

        for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
            NodeRef nr = arg.getNode(i);
            if (!arg.isRoot(nr)) {
                if (arg.getNodeHeight(nr) < time) {
                    Node left = (Node) arg.getParent(nr, 0);
                    Node right = (Node) arg.getParent(nr, 1);
                    if (arg.isBifurcation(nr)) {
                        assert left == right;
                        if (arg.getNodeHeight(left) > time) {
                            if (list != null)
                                list.add(nr);
                            count++;
                        }
                    } else {
                        if (arg.getNodeHeight(left) > time) {
                            if (list != null)
                                list.add(nr);
                            count++;
                        }
                        if (arg.getNodeHeight(right) > time) {
                            if (list != null)
                                list.add(nr);
                            count++;
                        }
                    }
                }
            } else {
                if (arg.getNodeHeight(nr) < time) {
                    if (list != null)
                        list.add(nr);
                    count++;
                }
            }
        }

        return count;
    }


    private int findPotentialNodesToRemove(ArrayList<NodeRef> list) {
        int count = 0;
        int n = arg.getNodeCount();


        for (int i = 0; i < n; i++) {
            Node node = (Node) arg.getNode(i);
            Node lp = node.leftParent;
            Node rp = node.rightParent;

            if (node.isReassortment() && (lp.bifurcation || rp.bifurcation)) {
                if (list != null)
                    list.add(node);
                count++;
            }
        }

        return count;
    }


    private double RemoveOperation() throws NoReassortmentEventException {

        double logHastings = 0;

        // 1. Draw reassortment node uniform randomly

        ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();
        int totalPotentials = findPotentialNodesToRemove(potentialNodes);

        if (totalPotentials == 0)
            throw new NoReassortmentEventException();

//		logHastings += Math.log((double)arg.getReassortmentNodeCount());
        logHastings += Math.log((double) totalPotentials);

//		double diff =(double)arg.getReassortmentNodeCount() - totalPotentials;
//
//		if(MathUtils.nextDouble() < diff/totalPotentials)
//			throw new NoReassortmentEventException();


        Node recNode = (Node) potentialNodes.get(MathUtils.nextInt(totalPotentials));

        double[] removePartitioningValues = recNode.partitioning.getParameterValues();

        double beforeReassortmentHeight = recNode.getHeight();
        double beforeBifurcationHeight = 0;
        double beforeTreeHeight = arg.getNodeHeight(arg.getRoot());

        arg.beginTreeEdit();

        boolean doneSomething = false;
        Node recParent = recNode.leftParent;
        Node recChild = recNode.leftChild;

        Node beforeReassortChild = recNode.leftChild;
        Node beforeBifurcationChild = recNode.leftParent;


        if (recNode.leftParent == recNode.rightParent) {
            if (!arg.isRoot(recNode.leftParent)) {
                beforeBifurcationHeight = recParent.getHeight();

                Node recGrandParent = recParent.leftParent;

                arg.doubleRemoveChild(recGrandParent, recParent);
                arg.doubleRemoveChild(recNode, recChild);
                if (recGrandParent.bifurcation)
                    arg.singleAddChild(recGrandParent, recChild);
                else
                    arg.doubleAddChild(recGrandParent, recChild);
                doneSomething = true;

                beforeBifurcationChild = beforeReassortChild;
            } else {
                //You should never go here.

                assert recChild.bifurcation;
                assert false;
            }
            logHastings += LOG_TWO;

        } else {

            Node recDeleteParent = null;
            Node recKeepParent = null;


            if (recNode.leftParent.bifurcation && recNode.rightParent.bifurcation) {
                if (MathUtils.nextBoolean()) {
                    recDeleteParent = recNode.rightParent;
                    recKeepParent = recNode.leftParent;
                } else {
                    recDeleteParent = recNode.leftParent;
                    recKeepParent = recNode.rightParent;
                }
                logHastings += LOG_TWO;
            } else if (recNode.rightParent.bifurcation) {
                recDeleteParent = recNode.rightParent;
                recKeepParent = recNode.leftParent;
            } else {
                recDeleteParent = recNode.leftParent;
                recKeepParent = recNode.rightParent;
            }

            beforeBifurcationChild = recDeleteParent.leftChild;
            if (beforeBifurcationChild == recNode) {
                beforeBifurcationChild = recDeleteParent.rightChild;
            }

            beforeBifurcationHeight = recDeleteParent.getHeight();


            if (arg.isRoot(recDeleteParent)) {

                Node oldRoot = (Node) arg.getOtherChild(recDeleteParent, recNode);
                Node oldRootLeft = oldRoot.leftChild;
                Node oldRootRight = oldRoot.rightChild;

                if (oldRoot == recKeepParent) {

                    arg.singleRemoveChild(recDeleteParent, recNode);
                    arg.singleRemoveChild(recDeleteParent, oldRoot);
                    arg.singleRemoveChild(oldRoot, oldRootLeft);
                    arg.singleRemoveChild(oldRoot, oldRootRight);

                    arg.singleAddChild(recDeleteParent, oldRootLeft);
                    arg.singleAddChild(recDeleteParent, oldRootRight);

                    arg.singleRemoveChild(recDeleteParent, recNode);
                    arg.doubleRemoveChild(recNode, recChild);

                    arg.singleAddChild(recDeleteParent, recChild);

                    recDeleteParent.setHeight(oldRoot.getHeight());

                    recDeleteParent = oldRoot;

                } else {
                    arg.singleRemoveChild(recDeleteParent, recNode);
                    arg.singleRemoveChild(recDeleteParent, oldRoot);
                    arg.singleRemoveChild(oldRoot, oldRootLeft);
                    arg.singleRemoveChild(oldRoot, oldRootRight);

                    arg.singleAddChild(recDeleteParent, oldRootLeft);
                    arg.singleAddChild(recDeleteParent, oldRootRight);

                    if (recKeepParent.bifurcation)
                        arg.singleRemoveChild(recKeepParent, recNode);
                    else
                        arg.doubleRemoveChild(recKeepParent, recNode);

                    arg.doubleRemoveChild(recNode, recChild);

                    if (recKeepParent.bifurcation)
                        arg.singleAddChild(recKeepParent, recChild);
                    else
                        arg.doubleAddChild(recKeepParent, recChild);

                    recDeleteParent.setHeight(oldRoot.getHeight());
                    recDeleteParent = oldRoot;

                }


            } else {
                Node recGrandParent = recDeleteParent.leftParent;

                Node otherChild = recDeleteParent.leftChild;
                if (otherChild == recNode)
                    otherChild = recDeleteParent.rightChild;

                if (recGrandParent.bifurcation)
                    arg.singleRemoveChild(recGrandParent, recDeleteParent);
                else
                    arg.doubleRemoveChild(recGrandParent, recDeleteParent);

                arg.singleRemoveChild(recDeleteParent, otherChild);
                if (recKeepParent.bifurcation)
                    arg.singleRemoveChild(recKeepParent, recNode);
                else
                    arg.doubleRemoveChild(recKeepParent, recNode);
                arg.doubleRemoveChild(recNode, recChild);
                if (otherChild != recChild) {
                    if (recGrandParent.bifurcation)
                        arg.singleAddChild(recGrandParent, otherChild);
                    else
                        arg.doubleAddChild(recGrandParent, otherChild);
                    if (recKeepParent.bifurcation)
                        arg.singleAddChild(recKeepParent, recChild);
                    else
                        arg.doubleAddChild(recKeepParent, recChild);
                } else {
                    if (recGrandParent.bifurcation)
                        arg.singleAddChildWithOneParent(recGrandParent, otherChild);
                    else
                        arg.doubleAddChildWithOneParent(recGrandParent, otherChild);
                    if (recKeepParent.bifurcation)
                        arg.singleAddChildWithOneParent(recKeepParent, recChild);
                    else
                        arg.doubleAddChildWithOneParent(recKeepParent, recChild);
                }
            }

            doneSomething = true;
            recParent = recDeleteParent;
        }

        if (relaxed) {

            double[] rateValues = recParent.rateParameter.getParameterValues();
            logHastings -= ratePrior.getAddHastingsRatio(rateValues);
        }


        if (doneSomething) {
            try {
                arg.contractARG(recParent, recNode,
                        internalNodeParameters, internalAndRootNodeParameters, nodeRates);

//                arg.contractARGWithRecombinant(recParent, recNode,
//                        internalNodeParameters, internalAndRootNodeParameters, nodeRates);
            } catch (Exception e) {
                System.err.println(e.getMessage());
                System.err.println(e);
                System.exit(-1);
            }
        }


        int max = Math.max(recParent.getNumber(), recNode.getNumber());
        int min = Math.min(recParent.getNumber(), recNode.getNumber());

        for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
            Node x = (Node) arg.getNode(i);
            if (x.getNumber() > max) {
                x.number--;
            }
            if (x.getNumber() > min) {
                x.number--;
            }
        }


        adjustRandomPartitioning();


        arg.pushTreeSizeDecreasedEvent();

        arg.endTreeEdit();
        try {
            arg.checkTreeIsValid();
        } catch (MutableTree.InvalidTreeException ite) {
            throw new RuntimeException(ite.toString() + "\n" + arg.toString()
                    + "\n" + TreeUtils.uniqueNewick(arg, arg.getRoot()));
        }

        assert nodeCheck() : arg.toARGSummary();

        //Do the backwards stuff now :(

        double afterTreeHeight = arg.getNodeHeight(arg.getRoot());

//		This is the ugly mixture proposal
//		double meanRoot = 4.0 / afterTreeHeight;
//		double case1 = 0.95;
//
//		if(beforeBifurcationHeight < afterTreeHeight){
//			logHastings -= 2.0*Math.log(afterTreeHeight) - Math.log(2.0*case1);
//		}else{
//			double additional = beforeBifurcationHeight - afterTreeHeight;
//			logHastings -= Math.log(afterTreeHeight) + additional*meanRoot -
//				Math.log((1-case1)*meanRoot);
//		}

        double theta = probBelowRoot / afterTreeHeight;

        logHastings -= theta * (beforeBifurcationHeight + beforeReassortmentHeight) - LOG_TWO
                - 2.0 * Math.log(theta) + Math.log(1 - Math.exp(-2.0 * afterTreeHeight * theta));


        logHastings -= Math.log((double) findPotentialAttachmentPoints(beforeBifurcationHeight, null)
                * findPotentialAttachmentPoints(beforeReassortmentHeight, null));


        logHastings -= getPartitionAddHastingsRatio(removePartitioningValues);

        logHastings -= LOG_TWO;

//        if (!beforeBifurcationChild.bifurcation &&
//                arg.getNodeHeight(beforeBifurcationChild.leftParent) > beforeBifurcationHeight &&
//                arg.getNodeHeight(beforeBifurcationChild.rightParent) > beforeBifurcationHeight) {
//            logHastings -= LOG_TWO;
//
//        }

//        if (!beforeReassortChild.bifurcation &&
//                arg.getNodeHeight(beforeReassortChild.leftParent) > beforeReassortmentHeight &&
//                arg.getNodeHeight(beforeReassortChild.rightParent) > beforeReassortmentHeight) {
//            logHastings -= LOG_TWO;
//        }
        assert nodeCheck();
        assert !Double.isNaN(logHastings) && !Double.isInfinite(logHastings);

        return logHastings;
    }

    private double getPartitionAddHastingsRatio(double[] values) {
        return -partLike.getLogLikelihood(values);
    }

    private void adjustRandomPartitioning() {
        if (tossSize < 1) {
            return;
        }

        if (arg.getReassortmentNodeCount() > 0) {
            int total = arg.getReassortmentNodeCount();
            Parameter xyz = arg.getPartitioningParameters().getParameter(MathUtils.nextInt(total));


            if (arg.isRecombinationPartitionType()) {
                adjustRecombinationPartition(xyz);
            } else {
                adjustReassortmentPartition(xyz);
            }
        }

    }

    private void adjustRecombinationPartition(Parameter part) {
        double[] values = part.getParameterValues();

        Logger.getLogger("dr.evomodel").severe("NOT IMPLENTED");
    }

    public static double arraySum(double[] n) {
        double b = 0;
        for (double a : n)
            b += a;
        return b;
    }

    private void adjustReassortmentPartition(Parameter part) {
        double[] values = part.getParameterValues();

        boolean stop = false;

        if (values.length == 2) {
            // Only option is to swap partitions
            double tmp = values[0];
            values[0] = values[1];
            values[1] = tmp;
        } else {

            while (!stop) {
                values = part.getParameterValues();

                ArrayList<Integer> list = new ArrayList<Integer>();

                while (list.size() < tossSize) {
                    int z = MathUtils.nextInt(values.length - 1) + 1;
                    if (!list.contains(z)) {
                        list.add(z);
                    }
                }


                for (int z : list) {
                    if (values[z] == 0) {
                        values[z] = 1;
                    } else {
                        values[z] = 0;
                    }
                }


                if (arraySum(values) > 0) {
                    stop = true;
                }
            }
        }

        for (int i = 0; i < values.length; i++) {


            part.setParameterValueQuietly(i, values[i]);
        }

        ARGPartitioningOperator.checkValidReassortmentPartition(part);


    }

    private void drawRandomPartitioning(Parameter partitioning) {
        double[] values = partLike.generatePartition();

        for (int i = 0; i < values.length; i++)
            partitioning.setParameterValueQuietly(i, values[i]);

    }

//    private int arraySum(int[] n) {
//        int b = 0;
//        for (int a : n)
//            b += a;
//        return b;
//    }
//
//
//
//    private void drawReassortmentPartitionAllFlip(Parameter partition) {
//        int numberOfPartitionsMinusOne = arg.getNumberOfPartitions() - 1;
//
//        int[] n = new int[numberOfPartitionsMinusOne];
//
//        while (arraySum(n) == 0) {
//            for (int i = 0; i < n.length; i++) {
//                if (MathUtils.nextBoolean()) {
//                    n[i] = 1;
//                } else {
//                    n[i] = 0;
//                }
//            }
//        }
//
//        partition.setParameterValueQuietly(0, 0);
//
//        for (int i = 0; i < numberOfPartitionsMinusOne; i++) {
//            partition.setParameterValueQuietly(i + 1, n[i]);
//        }
//
//
//
//        assert ARGPartitioningOperator.checkValidReassortmentPartition(partition);
//
//    }
//
//    private int nextFlipSize() {
//        int x = Poisson.nextPoisson(flipMean);
//        while (x > arg.getNumberOfPartitions() - 1 || x == 0) {
//            x = Poisson.nextPoisson(flipMean);
//        }
//
//        assert x != arg.getNumberOfPartitions();
//        assert x != 0;
//
//
//        return x;
//    }
//
//    private void drawReassortmentPartitionFixedFlip(Parameter partition) {
//        int numberOfPartitionsMinusOne = arg.getNumberOfPartitions() - 1;
//        int flipSize = nextFlipSize();
//
//        if (flipSize == numberOfPartitionsMinusOne) {
//            for (int i = 1; i <= flipSize; i++)
//                partition.setParameterValueQuietly(i, 1.0);
//            return;
//        }
//
//        double replaceValue = 1.0;
//        if (MathUtils.nextBoolean()) {
//            replaceValue = 0.0;
//
//            for (int i = 1; i <= numberOfPartitionsMinusOne; i++) {
//                partition.setParameterValueQuietly(i, 1.0);
//            }
//        }
//
//
//        ArrayList<Integer> a = new ArrayList<Integer>(flipSize);
//        while (a.size() < flipSize) {
//            int b = MathUtils.nextInt(numberOfPartitionsMinusOne) + 1;
//            if (!a.contains(b)) {
//                a.add(b);
//            }
//        }
//
//
//        for (int b : a) {
//            partition.setParameterValueQuietly(b, replaceValue);
//        }
//
//    }
//
//    private void drawRecombinationPartition(Parameter partition) {
//        int lengthMinusOne = arg.getNumberOfPartitions() - 1;
//
//        int cut = MathUtils.nextInt(lengthMinusOne);
//
//        int leftValue = 0;  //At one time, these values could switch.
//        int rightValue = 1;
//
//        for (int i = 0; i < cut + 1; i++)
//            partition.setParameterValueQuietly(i, leftValue);
//        for (int i = cut + 1; i < partition.getDimension(); i++)
//            partition.setParameterValueQuietly(i, rightValue);
//    }

    /* Draws a new partitioning.
      * With probability singlePartitionProbability, one bit is set;
      * otherwise, all bits are selected via a random permutation
      *
      */
    /* private double drawRandomReassortment(Parameter partitioning) {
         int len = arg.getNumberOfPartitions();
         double logq = 0;
         if (MathUtils.nextDouble() < singlePartitionProbability) {
             if (partitioning != null)
                 partitioning.setParameterValueQuietly(MathUtils.nextInt(len), 1.0);
             return Math.log(len);
         }
         int[] permutation = MathUtils.permuted(len);
         int cut = MathUtils.nextInt(len - 1);
         for (int i = 0; i < len; i++) {
             logq += Math.log(i + 1);
             if (i > cut && partitioning != null)
                 partitioning.setParameterValueQuietly(permutation[i], 1.0);
         }
         logq += Math.log(len - 1);
         return logq;
     } */

    /* Draws a new partitioning.
      * A break-pt is drawn uniformly
      *
      */
    /* private double drawRandomRecombination(Parameter partitioning) {
         int len = arg.getNumberOfPartitions();
         double logq = 0;
         double leftValue = MathUtils.nextInt(2);
         double rightValue = 1.0 - leftValue;
         logq += Math.log(2);
         if (partitioning != null) {
             int cut = MathUtils.nextInt(len - 1);
             for (int i = 0; i <= cut; i++)
                 partitioning.setParameterValueQuietly(i, leftValue);
             for (int i = cut + 1; i < len; i++)
                 partitioning.setParameterValueQuietly(i, rightValue);
         }
         logq += Math.log(len - 1);
         return logq;
     } */


    public boolean nodeCheck() {
        for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
            Node x = (Node) arg.getNode(i);

            if (x.leftParent != x.rightParent &&
                    x.leftChild != x.rightChild) {
                return false;
            }
            if (x.leftParent != null) {
                if (x.leftParent.leftChild.getNumber() != i &&
                        x.leftParent.rightChild.getNumber() != i)
                    return false;
            }
            if (x.rightParent != null) {
                if (x.rightParent.leftChild.getNumber() != i &&
                        x.rightParent.rightChild.getNumber() != i)
                    return false;
            }
            if (x.leftChild != null) {
                if (x.leftChild.leftParent.getNumber() != i &&
                        x.leftChild.rightParent.getNumber() != i)
                    return false;
            }
            if (x.rightChild != null) {
                if (x.rightChild.leftParent.getNumber() != i &&
                        x.rightChild.rightParent.getNumber() != i)
                    return false;
            }
        }

        return true;
    }


    public void sanityCheck() {
        int len = arg.getNodeCount();
        for (int i = 0; i < len; i++) {
            Node node = (Node) arg.getNode(i);
            if (node.bifurcation) {
                boolean equalChild = (node.leftChild == node.rightChild);
                if ((equalChild && node.leftChild != null)) {
                    if (!node.leftChild.bifurcation && ((node.leftChild).leftParent == node))
                        ;
                    else {
                        System.err.println("Node " + (i + 1) + " is insane.");
                        System.err.println(arg.toGraphString());
                        System.exit(-1);
                    }
                }
            } else {
                if ((node.leftChild != node.rightChild)) {
                    System.err.println("Node " + (i + 1) + " is insane.");
                    System.err.println(arg.toGraphString());
                    System.exit(-1);
                }
            }
            if (!node.isRoot()) {
                double d;
                d = node.getHeight();
            }
        }
    }


    ////
    ////Coercible MCMC Operator stuff
    ////


    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    public double getCoercableParameter() {
        return size;
    }

    public void setCoercableParameter(double value) {
        setSize(value);
    }

    public double getRawParameter() {
        return size;
    }

//	public int getMode() {
//		return mode;
//	}

    public double getTargetAcceptanceProbability() {
        return 0.5;
    }


    public String getPerformanceSuggestion() {
        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//		double targetProb = getTargetAcceptanceProbability();

//		double ws = OperatorUtils.optimizeWindowSize(getSize(), Double.MAX_VALUE, prob, targetProb);


        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting addProbability closer to 0.5";
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting addProbability value closer to 0.5";
        } else return "";
    }

    public String getOperatorName() {
        return ARG_EVENT_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return ARG_EVENT_OPERATOR;
        }

        //public String[] getParserNames(){
        //	return new String[]{getParserName(),"argEventOperator",};
        //}

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//			int mode = CoercableMCMCOperator.DEFAULT;
//			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
//				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
//					mode = CoercableMCMCOperator.COERCION_ON;
//				} else {
//					mode = CoercableMCMCOperator.COERCION_OFF;
//				}
//			}

            CoercionMode mode = CoercionMode.parseMode(xo);

            ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);

            CompoundParameter parameter1 = (CompoundParameter) xo.getChild(INTERNAL_NODES).getChild(0);
            CompoundParameter parameter2 = (CompoundParameter) xo.getChild(INTERNAL_AND_ROOT).getChild(0);
            CompoundParameter parameter3 = (CompoundParameter) xo.getChild(NODE_RATES).getChild(0);

            int weight = xo.getIntegerAttribute("weight");
            double size = xo.getDoubleAttribute(ADD_PROBABILITY);
            if (size > 0 && size < 1)
                size = Math.log(size / (1.0 - size));
            else
                throw new XMLParseException(ADD_PROBABILITY + " must be between 0 and 1");

            double belowRootProb = 0.9;
            if (xo.hasAttribute(BELOW_ROOT_PROBABILITY)) {
                belowRootProb = xo.getDoubleAttribute(BELOW_ROOT_PROBABILITY);
                if (belowRootProb >= 1 || belowRootProb <= 0) {
                    throw new XMLParseException(BELOW_ROOT_PROBABILITY + " must fall in (0,1)");
                }
            }

            ARGPartitionLikelihood partitionLike = (ARGPartitionLikelihood) xo.getChild(ARGPartitionLikelihood.class);

            int tossSize = 0;
            if (xo.hasAttribute(ARGPartitioningOperator.TOSS_SIZE)) {
                tossSize = xo.getIntegerAttribute(ARGPartitioningOperator.TOSS_SIZE);
                Logger.getLogger("dr.evomodel").info(ARG_EVENT_OPERATOR + " is joint with " + ARGPartitioningOperator.OPERATOR_NAME);
                if (tossSize <= 0 || tossSize >= treeModel.getNumberOfPartitions()) {
                    throw new XMLParseException("Toss size is incorrect");
                }
            }

            ARGRatePrior ratePrior = null;

            if (xo.hasAttribute(RELAXED)) {
                ratePrior = (ARGRatePrior) xo.getChild(ARGRatePrior.class);
            }

            return new ARGAddRemoveEventOperator(treeModel, weight, size,
                    mode, parameter1, parameter2, parameter3,
                    belowRootProb, partitionLike, ratePrior, tossSize);
        }

        public String getParserDescription() {
            return "An operator that slides a subarg.";
        }

        public Class getReturnType() {
            return ObsoleteARGAddRemoveEventOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(ADD_PROBABILITY, false,
                            "The probability that the operator adds a new"
                                    + " reassortment event"),
                    AttributeRule.newBooleanRule(JOINT_PARTITIONING, true),
                    AttributeRule.newIntegerRule(ARGPartitioningOperator.TOSS_SIZE, true),
                    AttributeRule.newDoubleRule(ADD_PROBABILITY, true),
                    AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                    AttributeRule.newIntegerRule(WEIGHT, false),
                    new ElementRule(ARGModel.class, false),
                    new ElementRule(ARGPartitionLikelihood.class, false),
                    new ElementRule(ARGRatePrior.class, true),
                    new ElementRule(INTERNAL_NODES,
                            new XMLSyntaxRule[]{
                                    new ElementRule(CompoundParameter.class)}),
                    new ElementRule(INTERNAL_AND_ROOT,
                            new XMLSyntaxRule[]{
                                    new ElementRule(CompoundParameter.class)}),
                    new ElementRule(NODE_RATES,
                            new XMLSyntaxRule[]{
                                    new ElementRule(CompoundParameter.class)}),
            };
        }
    };

    private class NoReassortmentEventException extends Exception {
        public NoReassortmentEventException(String message) {
            super(message);
        }

        public NoReassortmentEventException() {
            super("");
        }

        private static final long serialVersionUID = 1L;

    }

}

