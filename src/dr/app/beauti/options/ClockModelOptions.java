/*
 * ClockModelOptions.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.FixRateType;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.RelativeRatesType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import dr.evolution.util.Taxa;
import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.Vector;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class ClockModelOptions extends ModelOptions {

    private static final long serialVersionUID = 3544930558477534541L;
    // Instance variables
    private final BeautiOptions options;

//    private FixRateType rateOptionClockModel = FixRateType.RELATIVE_TO;
//    private double meanRelativeRate = 1.0;

//    public List<ClockModelGroup> clockModelGroupList = new ArrayList<ClockModelGroup>();

    public ClockModelOptions(BeautiOptions options) {
        this.options = options;

        initGlobalClockModelParaAndOpers();

//        fixRateOfFirstClockPartition();
    }

    private void initGlobalClockModelParaAndOpers() {

    }

    /**
     * return a list of parameters that are required
     */
    public void selectParameters() {
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
    }


    public String statusMessageClockModel(PartitionClockModel clockModel) {
        String t;
        if (clockModel.getRateTypeOption() == FixRateType.RELATIVE_TO) {
            if (clockModel.isEstimatedRate()) {
                t = "Estimate clock rate";
            } else {
                t = "Fix clock rate to " + clockModel.getRate();
            }

        } else {
            t = clockModel.getRateTypeOption().toString();
        }

        return t;
    }

    public double[] calculateInitialRootHeightAndRate(List<AbstractPartitionData> partitions) {
        double avgInitialRootHeight = 1;
        double avgInitialRate = 1;
        double avgMeanDistance = 1;

        if (partitions.size() > 0) {
            avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
        }

        if (options.getPartitionClockModels(partitions).size() > 0) {
            avgInitialRate = options.clockModelOptions.getSelectedRate(partitions); // all clock models
            PartitionClockModel clockModel = options.getPartitionClockModels(partitions).get(0);

            switch (clockModel.getRateTypeOption()) {
                case FIXED_MEAN:
                case RELATIVE_TO:
                    if (partitions.size() > 0) {
                        avgInitialRootHeight = avgMeanDistance / avgInitialRate;
                    }
                    break;

                case TIP_CALIBRATED:
                    avgInitialRootHeight = options.maximumTipHeight * 10.0;//TODO
                    avgInitialRate = avgMeanDistance / avgInitialRootHeight;//TODO
                    break;

                case NODE_CALIBRATED:
                    avgInitialRootHeight = getCalibrationEstimateOfRootTime(partitions);
                    if (avgInitialRootHeight < 0) avgInitialRootHeight = 1; // no leaf nodes
                    avgInitialRate = avgMeanDistance / avgInitialRootHeight;//TODO
                    break;

                case RATE_CALIBRATED:

                    break;

                default:
                    throw new IllegalArgumentException("Unknown fix rate type");
            }
        }
        avgInitialRootHeight = MathUtils.round(avgInitialRootHeight, 2);
        avgInitialRate = MathUtils.round(avgInitialRate, 2);

        return new double[]{avgInitialRootHeight, avgInitialRate};
    }

    public double getSelectedRate(List<AbstractPartitionData> partitions) {
        double selectedRate = 1;
        double avgInitialRootHeight;
        double avgMeanDistance = 1;

        if (partitions.size() > 0 && options.getPartitionClockModels(partitions).size() > 0) {
            //todo multi-group?
            PartitionClockModel clockModel = options.getPartitionClockModels(partitions).get(0);
            switch (clockModel.getRateTypeOption()) {
                case FIXED_MEAN:
                    selectedRate = clockModel.getRate();
                    break;

                case RELATIVE_TO:
                    List<PartitionClockModel> models = options.getPartitionClockModels(partitions);
                    // fix ?th partition
                    if (models.size() == 1) {
                        selectedRate = models.get(0).getRate();
                    } else {
                        selectedRate = getAverageRate(models);
                    }
                    break;

                case TIP_CALIBRATED:
                    if (partitions.size() > 0) {
                        avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
                    }
                    avgInitialRootHeight = options.maximumTipHeight * 10.0;//TODO
                    selectedRate = avgMeanDistance / avgInitialRootHeight;//TODO
                    break;

                case NODE_CALIBRATED:
                    if (partitions.size() > 0) {
                        avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
                    }
                    avgInitialRootHeight = getCalibrationEstimateOfRootTime(partitions);
                    if (avgInitialRootHeight < 0) avgInitialRootHeight = 1; // no leaf nodes
                    selectedRate = avgMeanDistance / avgInitialRootHeight;//TODO
                    break;

                case RATE_CALIBRATED:
                    //TODO
                    break;

                default:
                    throw new IllegalArgumentException("Unknown fix rate type");
            }
        }
        return selectedRate;
    }

    private double getCalibrationEstimateOfRootTime(List<AbstractPartitionData> partitions) {

        // TODO - shouldn't this method be in the PartitionTreeModel??

        List<Taxa> taxonSets = options.taxonSets;
        if (taxonSets != null && taxonSets.size() > 0) { // tmrca statistic

            // estimated root times based on each of the taxon sets
            double[] rootTimes = new double[taxonSets.size()];

            for (int i = 0; i < taxonSets.size(); i++) {

                Taxa taxa = taxonSets.get(i);

                Parameter tmrcaStatistic = options.getStatistic(taxa);

                double taxonSetCalibrationTime = tmrcaStatistic.getPriorExpectationMean();

                // the calibration distance is the patristic genetic distance back to the common ancestor of
                // the set of taxa.
                double calibrationDistance = 0;

                // the root distance is the patristic genetic distance back to the root of the tree.
                double rootDistance = 0;

                int siteCount = 0;

                for (AbstractPartitionData partition : partitions) {
                    if (partition.getDistances() != null) {   // ignore partitions that don't have distances
                        Tree tree = new UPGMATree(partition.getDistances());

                        Set<String> leafNodes = Taxa.Utils.getTaxonListIdSet(taxa);

                        if (leafNodes.size() < 1) {
                            return -1;
                        }

                        NodeRef node = Tree.Utils.getCommonAncestorNode(tree, leafNodes);

                        if (node == null) {
                            throw new IllegalArgumentException("Can't find MRCA node for taxon set, " + taxa.getId() + ", in partition: " + partition.getName());
                        }

                        calibrationDistance += tree.getNodeHeight(node);
                        rootDistance += tree.getNodeHeight(tree.getRoot());

                        siteCount += partition.getSiteCount();
                    }
                }

                rootDistance /= partitions.size();
                calibrationDistance /= partitions.size();

                if (calibrationDistance == 0.0) {
                    calibrationDistance = 0.25 / siteCount;
                }

                if (rootDistance == 0) {
                    rootDistance = 0.25 / siteCount;
                }

                rootTimes[i] += (rootDistance / calibrationDistance) * taxonSetCalibrationTime;
            }

            // return the mean estimate of the root time for this set of partitions
            return DiscreteStatistics.mean(rootTimes);

        } else { // prior on treeModel.rootHight
            double avgInitialRootHeight = 0;
            double count = 0;
            for (PartitionTreeModel tree : options.getPartitionTreeModels(partitions)) {
                avgInitialRootHeight = avgInitialRootHeight + tree.getInitialRootHeight();
                count = count + 1;
            }
            if (count != 0) avgInitialRootHeight = avgInitialRootHeight / count;
            return avgInitialRootHeight;
        }

    }

    public double getAverageRate(List<PartitionClockModel> models) { //TODO average per tree, but how to control the estimate clock => tree?
        double averageRate = 0;
        double count = 0;

        for (PartitionClockModel model : models) {
            if (!model.isEstimatedRate()) {
                averageRate = averageRate + model.getRate();
                count = count + 1;
            }
        }

        if (count > 0) {
            averageRate = averageRate / count;
        } else {
            averageRate = 1; //TODO how to calculate rate when estimate all
        }

        return averageRate;
    }

    public boolean isTipCalibrated() {
        return options.maximumTipHeight > 0;
    }

    public boolean isRateCalibrated() {
        return false;//TODO
    }

    public int[] getPartitionClockWeights() {
        int[] weights = new int[options.getPartitionClockModels().size()]; // use List?

        int k = 0;
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            for (AbstractPartitionData partition : options.getDataPartitions(model)) {
                int n = partition.getSiteCount();
                weights[k] += n;
            }
            k += 1;
        }

        assert (k == weights.length);

        return weights;
    }

}
