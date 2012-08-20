/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
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

//        createParameter("allClockRates", "All the relative rates regarding clock models");
//
//        createOperator("deltaAllClockRates", RelativeRatesType.CLOCK_RELATIVE_RATES.toString(),
//                "Delta exchange operator for all the relative rates regarding clock models", "allClockRates",
//                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);
//
//        // only available for *BEAST and EBSP
//        createUpDownAllOperator("upDownAllRatesHeights", "Up down all rates and heights",
//                "Scales all rates inversely to node heights of the tree",
//                demoTuning, branchWeights);

    }

    /**
     * return a list of parameters that are required
     */
    public void selectParameters() {
        for (ClockModelGroup clockModelGroup : getClockModelGroups()) {
            createParameter(clockModelGroup.getName(), // used in BeastGenerator (Branch Rates Model) part
                    "All relative rates regarding clock models in group " + clockModelGroup.getName());
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        for (ClockModelGroup clockModelGroup : getClockModelGroups()) {
            if (clockModelGroup.getRateTypeOption() == FixRateType.FIX_MEAN) {
                createOperator("delta_" + clockModelGroup.getName(),
                        RelativeRatesType.CLOCK_RELATIVE_RATES.toString() + " in " + clockModelGroup.getName(),
                        "Delta exchange operator for all relative rates regarding clock models",
                        clockModelGroup.getName(), OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);

                Operator deltaOperator = getOperator("delta_" + clockModelGroup.getName());

                // update delta clock operator weight
                deltaOperator.weight = options.getPartitionClockModels(clockModelGroup).size();

                ops.add(deltaOperator);
            }

            //up down all rates and trees operator only available for *BEAST and EBSP
            if (clockModelGroup.getRateTypeOption() == FixRateType.RELATIVE_TO && //TODO what about Calibration?
                    (options.useStarBEAST || options.isEBSPSharingSamePrior())) {
                // only available for *BEAST and EBSP
                createUpDownAllOperator("upDownAllRatesHeights_" + clockModelGroup.getName(),
                        "Up down all rates and heights in " + clockModelGroup.getName(),
                        "Scales all rates inversely to node heights of the tree",
                        demoTuning, branchWeights);
                Operator op = getOperator("upDownAllRatesHeights_" + clockModelGroup.getName());
                op.setClockModelGroup(clockModelGroup);

                ops.add(op);
            }
        }
    }

    //+++++++++++++++++++++++ Clock Model Group ++++++++++++++++++++++++++++++++
    public void initClockModelGroup() { // only used in BeautiImporter
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            addClockModelGroup(model);
        }

        for (ClockModelGroup clockModelGroup : getClockModelGroups()) {
            if (clockModelGroup.contain(Microsatellite.INSTANCE, options)) {
                if (options.getPartitionClockModels(clockModelGroup).size() == 1) {
                    fixRateOfFirstClockPartition(clockModelGroup);
                    options.getPartitionClockModels(clockModelGroup).get(0).setEstimatedRate(true);
                } else {
                    fixMeanRate(clockModelGroup);
                }
            } else if (!(clockModelGroup.getRateTypeOption() == FixRateType.TIP_CALIBRATED
                    || clockModelGroup.getRateTypeOption() == FixRateType.NODE_CALIBRATED
                    || clockModelGroup.getRateTypeOption() == FixRateType.RATE_CALIBRATED)) {
                //TODO correct?
                fixRateOfFirstClockPartition(clockModelGroup);
            }
        }
    }

    public void addClockModelGroup(PartitionClockModel model) {
        if (model.getClockModelGroup() == null) {
            String groupName = model.getDataType().getDescription().toLowerCase() + "_group";
            List<ClockModelGroup> groupsList = getClockModelGroups();
            ClockModelGroup clockModelGroup;
            if (containsGroup(groupName, groupsList)) {
                clockModelGroup = getGroup(groupName, groupsList);
            } else {
                clockModelGroup = new ClockModelGroup(groupName);
            }
            model.setClockModelGroup(clockModelGroup);
        }
    }

    public List<ClockModelGroup> getClockModelGroups(DataType dataType) {
        List<ClockModelGroup> activeClockModelGroups = new ArrayList<ClockModelGroup>();
        for (PartitionClockModel model : options.getPartitionClockModels(dataType)) {
            ClockModelGroup group = model.getClockModelGroup();
            if (group != null && (!activeClockModelGroups.contains(group))) {
                activeClockModelGroups.add(group);
            }
        }
        return activeClockModelGroups;
    }

    public List<ClockModelGroup> getClockModelGroups(List<? extends AbstractPartitionData> givenDataPartitions) {
        List<ClockModelGroup> activeClockModelGroups = new ArrayList<ClockModelGroup>();
        for (PartitionClockModel model : options.getPartitionClockModels(givenDataPartitions)) {
            ClockModelGroup group = model.getClockModelGroup();
            if (group != null && (!activeClockModelGroups.contains(group))) {
                activeClockModelGroups.add(group);
            }
        }
        return activeClockModelGroups;
    }

    public List<ClockModelGroup> getClockModelGroups() {
        return getClockModelGroups(options.dataPartitions);
    }

    public Vector<String> getClockModelGroupNames(List<ClockModelGroup> group) {
        Vector<String> activeClockModelGroups = new Vector<String>();
        for (ClockModelGroup clockModelGroup : group) {
            String name = clockModelGroup.getName();
            if (name != null && (!activeClockModelGroups.contains(name))) {
                activeClockModelGroups.add(name);
            }
        }
        return activeClockModelGroups;
    }

    public boolean containsGroup(String groupName, List<ClockModelGroup> groupsList) {
        for (ClockModelGroup clockModelGroup : groupsList) {
            if (clockModelGroup.getName().equalsIgnoreCase(groupName)) return true;
        }
        return false;
    }

    public ClockModelGroup getGroup(String groupName, List<ClockModelGroup> groupsList) {
        for (ClockModelGroup clockModelGroup : groupsList) {
            if (clockModelGroup.getName().equalsIgnoreCase(groupName))
                return clockModelGroup;
        }
        return null;
    }

    public void fixRateOfFirstClockPartition(ClockModelGroup group) {
        group.setRateTypeOption(FixRateType.RELATIVE_TO);
        // fix rate of 1st partition
        int i = 0;
        for (PartitionClockModel model : options.getPartitionClockModels(group)) {
            if (i < 1) {
                model.setEstimatedRate(false);
            } else {
                model.setEstimatedRate(true);
            }
            i = i + 1;
        }
    }

    public void fixMeanRate(ClockModelGroup group) {
        group.setRateTypeOption(FixRateType.FIX_MEAN);

        for (PartitionClockModel model : options.getPartitionClockModels(group)) {
            model.setEstimatedRate(true); // all set to NOT fixed, because detla exchange
            model.setRate(group.getFixMeanRate(), false);
        }
    }

    public void tipTimeCalibration(ClockModelGroup group) {
        group.setRateTypeOption(FixRateType.TIP_CALIBRATED);

        for (PartitionClockModel model : options.getPartitionClockModels(group)) {
            model.setEstimatedRate(true);
        }
    }


    public void nodeCalibration(ClockModelGroup group) {
        group.setRateTypeOption(FixRateType.NODE_CALIBRATED);

        for (PartitionClockModel model : options.getPartitionClockModels(group)) {
            model.setEstimatedRate(true);
        }
    }


    public void rateCalibration(ClockModelGroup group) {
        group.setRateTypeOption(FixRateType.RATE_CALIBRATED);

        for (PartitionClockModel model : options.getPartitionClockModels(group)) {
            model.setEstimatedRate(true);
        }
    }

    public String statusMessageClockModel(ClockModelGroup group) {
        String t;
        if (group.getRateTypeOption() == FixRateType.RELATIVE_TO) {
            if (options.getPartitionClockModels(group).size() == 1) { // single partition clock
                if (options.getPartitionClockModels(group).get(0).isEstimatedRate()) {
                    t = "Estimate clock rate";
                } else {
                    t = "Fix clock rate to " + options.getPartitionClockModels(group).get(0).getRate();
                }

            } else {
                // todo is the following code excuted?
                t = group.getRateTypeOption().toString() + " ";
                int c = 0;
                for (PartitionClockModel model : options.getPartitionClockModels(group)) {
                    if (!model.isEstimatedRate()) {
                        if (c > 0) t = t + ", ";
                        c = c + 1;
                        t = t + model.getName();
                    }
                }

                if (c == 0) t = "Estimate all clock rates";
                if (c == options.getPartitionClockModels(group).size()) t = "Fix all clock rates";
            }

        } else {
            t = group.getRateTypeOption().toString();
        }

        return t + " in " + group.getName();
    }

    public String statusMessageClockModel() {
        String t = "";
        for (ClockModelGroup clockModelGroup : getClockModelGroups()) {
            t += statusMessageClockModel(clockModelGroup) + "; ";
        }
        return t;
    }


    /////////////////////////////////////////////////////////////
//    public FixRateType getRateOptionClockModel() {
//        return rateOptionClockModel;
//    }

//	public void setRateOptionClockModel(FixRateType rateOptionClockModel) {
//		this.rateOptionClockModel = rateOptionClockModel;
//	}

//    public void setMeanRelativeRate(double meanRelativeRate) {
//        this.meanRelativeRate = meanRelativeRate;
//    }

//    public double calculateAvgBranchLength(List<AbstractPartitionData> partitions) { // todo
//        double avgBranchLength = 1;
//
//        for (PartitionTreeModel tree : options.getPartitionTreeModels(partitions)) {
//        }
//        return MathUtils.round(avgBranchLength, 2);
//    }


    public double[] calculateInitialRootHeightAndRate(List<AbstractPartitionData> partitions) {
        double avgInitialRootHeight = 1;
        double avgInitialRate = 1;
        double avgMeanDistance = 1;

//        List<AbstractPartitionData> partitions = options.getDataPartitions(clockModelGroup);

        if (partitions.size() > 0) {
            avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
        }

        if (options.getPartitionClockModels(partitions).size() > 0) {
            avgInitialRate = options.clockModelOptions.getSelectedRate(partitions); // all clock models
            //todo multi-group?
            ClockModelGroup clockModelGroup = options.getPartitionClockModels(partitions).get(0).getClockModelGroup();

            switch (clockModelGroup.getRateTypeOption()) {
                case FIX_MEAN:
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
        // calibration: all isEstimatedRate = true

//        List<AbstractPartitionData> partitions = options.getDataPartitions(clockModelGroup);

        if (partitions.size() > 0 && options.getPartitionClockModels(partitions).size() > 0) {
            //todo multi-group?
            ClockModelGroup clockModelGroup = options.getPartitionClockModels(partitions).get(0).getClockModelGroup();
            switch (clockModelGroup.getRateTypeOption()) {
                case FIX_MEAN:
                    selectedRate = clockModelGroup.getFixMeanRate();
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

//    private List<AbstractPartitionData> getAllPartitionDataGivenClockModels(List<PartitionClockModel> models) {
//
//        List<AbstractPartitionData> allData = new ArrayList<AbstractPartitionData>();
//
//        for (PartitionClockModel model : models) {
//            for (AbstractPartitionData partition : model.getDataPartitions()) {
//                if (partition != null && (!allData.contains(partition))) {
//                    allData.add(partition);
//                }
//            }
//        }
//
//        return allData;
//    }

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
                    Tree tree = new UPGMATree(partition.getDistances());

                    Set<String> leafNodes = Taxa.Utils.getTaxonListIdSet(taxa);

                    if (leafNodes.size() < 1) {
                        return -1;
                    }

                    NodeRef node = Tree.Utils.getCommonAncestorNode(tree, leafNodes);

                    calibrationDistance += tree.getNodeHeight(node);
                    rootDistance += tree.getNodeHeight(tree.getRoot());

                    siteCount += partition.getSiteCount();
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

    // FixRateType.FIX_MEAN
//    public double getMeanRelativeRate() {
//        return meanRelativeRate;
//    }

    // FixRateType.ESTIMATE

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

    // Calibration Series Data
    public double getAverageRateForCalibrationSeriesData() {
        //TODO
        return (double) 0;
    }

    // Calibration TMRCA
    public double getAverageRateForCalibrationTMRCA() {
        //TODO
        return (double) 0;
    }

    public boolean isTipCalibrated() {
        return options.maximumTipHeight > 0;
    }

    public boolean isRateCalibrated() {
        return false;//TODO
    }

    public int[] getPartitionClockWeights(ClockModelGroup group) {
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

//    public void fixRateOfFirstClockPartition() {
//        this.rateOptionClockModel = FixRateType.RELATIVE_TO;
//        // fix rate of 1st partition
//        int i = 0;
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            if (i < 1) {
//                model.setEstimatedRate(false);
//            } else {
//                model.setEstimatedRate(true);
//            }
//            i = i + 1;
//        }
//    }
//
//    public void fixMeanRate() {
//        this.rateOptionClockModel = FixRateType.FIX_MEAN;
//
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            model.setEstimatedRate(true); // all set to NOT fixed, because detla exchange
//        }
//    }
//
//    public void tipTimeCalibration() {
//        this.rateOptionClockModel = FixRateType.TIP_CALIBRATED;
//
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            model.setEstimatedRate(true);
//        }
//    }
//
//
//    public void nodeCalibration() {
//        this.rateOptionClockModel = FixRateType.NODE_CALIBRATED;
//
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            model.setEstimatedRate(true);
//        }
//    }
//
//
//    public void rateCalibration() {
//        this.rateOptionClockModel = FixRateType.RATE_CALIBRATED;
//
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            model.setEstimatedRate(true);
//        }
//    }

//    public String statusMessageClockModel() {
//        if (rateOptionClockModel == FixRateType.RELATIVE_TO) {
//            if (options.getPartitionClockModels().size() == 1) { // single partition clock
//                if (options.getPartitionClockModels().get(0).isEstimatedRate()) {
//                    return "Estimate clock rate";
//                } else {
//                    return "Fix clock rate to " + options.getPartitionClockModels().get(0).getRate();
//                }
//
//            } else {
//                String t = rateOptionClockModel.toString() + " ";
//                int c = 0;
//                for (PartitionClockModel model : options.getPartitionClockModels()) {
//                    if (!model.isEstimatedRate()) {
//                        if (c > 0) t = t + ", ";
//                        c = c + 1;
//                        t = t + model.getName();
//                    }
//                }
//
//                if (c == 0) t = "Estimate all clock rates";
//                if (c == options.getPartitionClockModels().size()) t = "Fix all clock rates";
//
//                return t;
//            }
//
//        } else {
//            return rateOptionClockModel.toString();
//        }
//    }

    //+++++++++++++++++++++++ Validation ++++++++++++++++++++++++++++++++
    // true => valid, false => warning message
//    public boolean validateFixMeanRate(boolean fixedMeanRateCheck) {
//        return !(fixedMeanRateCheck && options.getPartitionClockModels().size() < 2);
//    }

//    public boolean validateRelativeTo() {
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            if (!model.isEstimatedRate()) { // fixed
//                return true;
//            }
//        }
//        return false;
//    }

}
