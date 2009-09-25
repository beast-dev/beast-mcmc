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

import dr.app.beauti.enumTypes.FixRateType;
import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorType;
import dr.app.beauti.enumTypes.RelativeRatesType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import dr.evolution.util.Taxa;
import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;

import java.util.List;

import javax.swing.JCheckBox;



/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class ClockModelOptions extends ModelOptions {

    // Instance variables
    private final BeautiOptions options;

    private FixRateType rateOptionClockModel = FixRateType.RELATIVE_TO;
    private double meanRelativeRate = 1.0;

    public ClockModelOptions(BeautiOptions options) {
        this.options = options;

        initGlobalClockModelParaAndOpers();

        fixRateOfFirstClockPartition();
    }

    private void initGlobalClockModelParaAndOpers() {

        createParameter("allClockRates", "All the relative rates regarding clock models");

        createOperator("deltaAllClockRates", RelativeRatesType.CLOCK_RELATIVE_RATES.toString(),
                "Delta exchange operator for all the relative rates regarding clock models", "allClockRates",
                OperatorType.DELTA_EXCHANGE, 0.75, rateWeights);

        // only available for *BEAST and EBSP
        createUpDownAllOperator("upDownAllRatesHeights", "Up down all rates and heights", "Scales all rates inversely to node heights of the tree",
                demoTuning, branchWeights);

    }

    /**
     * return a list of parameters that are required
     *
     */
    public void selectParameters() {

    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (rateOptionClockModel == FixRateType.FIX_MEAN) {
            Operator deltaOperator = getOperator("deltaAllClockRates");

            // update delta clock operator weight
            deltaOperator.weight = options.getPartitionClockModels().size();

            ops.add(deltaOperator);
        }

        //up down all rates and trees operator only available for *BEAST and EBSP
        if (rateOptionClockModel == FixRateType.RELATIVE_TO && //TODO what about Calibration? 
                (options.starBEASTOptions.isSpeciesAnalysis() || options.isEBSPSharingSamePrior())) {
            ops.add(getOperator("upDownAllRatesHeights"));
        }

    }


    /////////////////////////////////////////////////////////////
    public FixRateType getRateOptionClockModel() {
        return rateOptionClockModel;
    }

//	public void setRateOptionClockModel(FixRateType rateOptionClockModel) {
//		this.rateOptionClockModel = rateOptionClockModel;
//	}

    public void setMeanRelativeRate(double meanRelativeRate) {
        this.meanRelativeRate = meanRelativeRate;
    }

    public double[] calculateInitialRootHeightAndRate(List<PartitionData> partitions) {
        double avgInitialRootHeight = 1;
        double avgInitialRate = 1;
        double avgMeanDistance = 1;

        if (options.hasData()) {
            avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
        }

        if (options.getPartitionClockModels().size() > 0) {
            avgInitialRate = options.clockModelOptions.getSelectedRate(partitions); // all clock models
        }

        switch (options.clockModelOptions.getRateOptionClockModel()) {
            case FIX_MEAN:
            case RELATIVE_TO:
                if (options.hasData()) {
                    avgInitialRootHeight = avgMeanDistance / avgInitialRate;
                }
                break;

            case TIP_CALIBRATED:
                avgInitialRootHeight = options.maximumTipHeight * 10.0;//TODO
                avgInitialRate = avgMeanDistance / avgInitialRootHeight;//TODO
                break;

            case NODE_CALIBRATED:
            	avgInitialRootHeight = getCalibrationEstimateOfRootTime(partitions);
                avgInitialRate = avgMeanDistance / avgInitialRootHeight;//TODO
                break;

            case RATE_CALIBRATED:

                break;

            default:
                throw new IllegalArgumentException("Unknown fix rate type");
        }

        avgInitialRootHeight = MathUtils.round(avgInitialRootHeight, 2);
        avgInitialRate = MathUtils.round(avgInitialRate, 2);

        return new double[]{avgInitialRootHeight, avgInitialRate};
    }

    public double getSelectedRate(List<PartitionData> partitions) {
        double selectedRate = 1;
        double avgInitialRootHeight;
        double avgMeanDistance = 1;
        // calibration: all isEstimatedRate = true
        switch (options.clockModelOptions.getRateOptionClockModel()) {
            case FIX_MEAN:
                selectedRate = meanRelativeRate;
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
                if (options.hasData()) {
                    avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
                }
                avgInitialRootHeight = options.maximumTipHeight * 10.0;//TODO
                selectedRate = avgMeanDistance / avgInitialRootHeight;//TODO
                break;

            case NODE_CALIBRATED:
                if (options.hasData()) {
                    avgMeanDistance = options.getAveWeightedMeanDistance(partitions);
                }
                avgInitialRootHeight = getCalibrationEstimateOfRootTime(partitions);
                selectedRate = avgMeanDistance / avgInitialRootHeight;//TODO
                break;

            case RATE_CALIBRATED:
              //TODO
                break;

            default:
                throw new IllegalArgumentException("Unknown fix rate type");
        }
        
        return selectedRate;
    }
    
//    private List<PartitionData> getAllPartitionDataGivenClockModels(List<PartitionClockModel> models) {
//
//        List<PartitionData> allData = new ArrayList<PartitionData>();
//
//        for (PartitionClockModel model : models) {
//            for (PartitionData partition : model.getAllPartitionData()) {
//                if (partition != null && (!allData.contains(partition))) {
//                    allData.add(partition);
//                }
//            }
//        }
//
//        return allData;
//    }

    private double getCalibrationEstimateOfRootTime(List<PartitionData> partitions) {

        // TODO - shouldn't this method be in the PartitionTreeModel??

        List<Taxa> taxonSets = options.taxonSets;
        if (taxonSets != null && taxonSets.size() > 0) { // tmrca statistic 

            // estimated root times based on each of the taxon sets
            double[] rootTimes = new double[taxonSets.size()];

            for (int i = 0; i < taxonSets.size(); i++) {

                Taxa taxa = taxonSets.get(i);

                Parameter tmrcaStatistic = options.getStatistic(taxa);

                double taxonSetCalibrationTime = tmrcaStatistic.getPriorExpectation();

                // the calibration distance is the patristic genetic distance back to the common ancestor of
                // the set of taxa.
                double calibrationDistance = 0;

                // the root distance is the patristic genetic distance back to the root of the tree.
                double rootDistance = 0;

                int siteCount = 0;

                for (PartitionData partition : partitions) {
                    Tree tree = new UPGMATree(partition.distances);

                    NodeRef node = Tree.Utils.getCommonAncestorNode(tree, Taxa.Utils.getTaxonListIdSet(taxa));

                    calibrationDistance += tree.getNodeHeight(node);
                    rootDistance += tree.getNodeHeight(tree.getRoot());

                    siteCount += partition.getAlignment().getSiteCount();
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
    public double getMeanRelativeRate() {
        return meanRelativeRate;
    }

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

    public boolean isNodeCalibrated(Parameter para) {
        return (para.taxa != null && hasProperPriorOn(para)) // param.taxa != null is TMRCA
                || (para.getBaseName().endsWith("treeModel.rootHeight") && hasProperPriorOn(para));
    }
    
    private boolean hasProperPriorOn(Parameter para) {
        return para.priorType == PriorType.LOGNORMAL_PRIOR
                || para.priorType == PriorType.NORMAL_PRIOR
                || para.priorType == PriorType.LAPLACE_PRIOR
                || para.priorType == PriorType.TRUNC_NORMAL_PRIOR
                || (para.priorType == PriorType.GAMMA_PRIOR && para.shape > 1)
                || (para.priorType == PriorType.GAMMA_PRIOR && para.offset > 0)
                || (para.priorType == PriorType.UNIFORM_PRIOR && para.lower > 0 && para.upper < Double.POSITIVE_INFINITY)
                || (para.priorType == PriorType.EXPONENTIAL_PRIOR && para.offset > 0);
    }


    public boolean isRateCalibrated() {
        return false;//TODO
    }

    public int[] getPartitionClockWeights() {
        int[] weights = new int[options.getPartitionClockModels().size()]; // use List?

        int k = 0;
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            for (PartitionData partition : model.getAllPartitionData()) {
                int n = partition.getSiteCount();
                weights[k] += n;
            }
            k += 1;
        }

        assert (k == weights.length);

        return weights;
    }

    public void fixRateOfFirstClockPartition() {
        this.rateOptionClockModel = FixRateType.RELATIVE_TO;
        // fix rate of 1st partition
        int i = 0;
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            if (i < 1) {
                model.setEstimatedRate(false);
            } else {
                model.setEstimatedRate(true);
            }
            i = i + 1;
        }
    }

    public void fixMeanRate() {
        this.rateOptionClockModel = FixRateType.FIX_MEAN;

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            model.setEstimatedRate(true); // all set to NOT fixed, because detla exchange
        }
    }

    public void tipTimeCalibration() {
        this.rateOptionClockModel = FixRateType.TIP_CALIBRATED;

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            model.setEstimatedRate(true);
        }
    }


    public void nodeCalibration() {
        this.rateOptionClockModel = FixRateType.NODE_CALIBRATED;

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            model.setEstimatedRate(true);
        }
    }


    public void rateCalibration() {
        this.rateOptionClockModel = FixRateType.RATE_CALIBRATED;

        for (PartitionClockModel model : options.getPartitionClockModels()) {
            model.setEstimatedRate(true);
        }
    }
    
    public String statusMessageClockModel() {
        if (rateOptionClockModel == FixRateType.RELATIVE_TO) {
            if (options.getPartitionClockModels().size() == 1) { // single partition clock
                if (options.getPartitionClockModels().get(0).isEstimatedRate()) {
                    return "Estimate clock rate";
                } else {
                    return "Fix clock rate to " + options.getPartitionClockModels().get(0).getRate();
                }
                
            } else {
                String t = rateOptionClockModel.toString() + " ";
                int c = 0;
                for (PartitionClockModel model : options.getPartitionClockModels()) {
                    if (!model.isEstimatedRate()) {
                        if (c > 0) t = t + ", ";
                        c = c + 1;
                        t = t + model.getName();                    
                    }
                }
                
                if (c == 0) t = "Estimate all clock rates";
                if (c == options.getPartitionClockModels().size()) t = "Fix all clock rates";
                
                return t;                
            }
            
        } else {
            return rateOptionClockModel.toString();
        }
    }

    //+++++++++++++++++++++++ Validation ++++++++++++++++++++++++++++++++
    // true => valid, false => warning message 
    public boolean validateFixMeanRate(JCheckBox fixedMeanRateCheck) {
        return !(fixedMeanRateCheck.isSelected() && options.getPartitionClockModels().size() < 2);
    }
    
//    public boolean validateRelativeTo() {
//        for (PartitionClockModel model : options.getPartitionClockModels()) {
//            if (!model.isEstimatedRate()) { // fixed
//                return true;
//            }
//        }        
//        return false;
//    }

}
