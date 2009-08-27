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
import javax.swing.JOptionPane;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class ClockModelOptions extends ModelOptions {

    // Instance variables
    private final BeautiOptions options;

    private FixRateType rateOptionClockModel = FixRateType.RElATIVE_TO;
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
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {

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
        if (rateOptionClockModel == FixRateType.RElATIVE_TO && //TODO what about Calibration? 
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
            avgInitialRate = options.clockModelOptions.getSelectedRate(options.getPartitionClockModels(partitions)); // all clock models
        }

        switch (options.clockModelOptions.getRateOptionClockModel()) {
            case FIX_MEAN:
            case RElATIVE_TO:
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

    public double getSelectedRate(List<PartitionClockModel> models) {
        double selectedRate = 1;
        if (rateOptionClockModel == FixRateType.FIX_MEAN) {
            selectedRate = meanRelativeRate;

        } else if (rateOptionClockModel == FixRateType.RElATIVE_TO) {
            // fix ?th partition
            if (models.size() == 1) {
                selectedRate = models.get(0).getRate();
            } else {
                selectedRate = getAverageRate(models);
            }

        } else {
            // calibration: all isEstimatedRate = true
            //TODO calibration
        }
        return selectedRate;
    }

    private double getCalibrationEstimateOfRootTime(List<PartitionData> partitions) {

        // TODO - shouldn't this method be in the PartitionTreeModel??

        List<Taxa> taxonSets = options.taxonSets;
        if (taxonSets != null && taxonSets.size() > 0) {

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
        }

        return 0.0;
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
        }

        return averageRate;
    }

    // Calibration Series Data
    public double getAverageRateForCalibrationSeriesData() {
        double averageRate = 0;
        //TODO
        return averageRate;
    }

    // Calibration TMRCA
    public double getAverageRateForCalibrationTMRCA() {
        double averageRate = 0;
        //TODO
        return averageRate;
    }

    public boolean isTipCalibrated() {
        return options.maximumTipHeight > 0;
    }

    public boolean isNodeCalibrated(Parameter para) {
        if ((para.taxa != null && hasProperPriorOn(para)) // param.taxa != null is TMRCA
                || (para.getBaseName().endsWith("treeModel.rootHeight") && hasProperPriorOn(para))) {
            return true;
        } else {            
            return false;
        }
    }
    
    private boolean hasProperPriorOn(Parameter para) {                
        if (para.priorType == PriorType.LOGNORMAL_PRIOR
                || para.priorType == PriorType.NORMAL_PRIOR
                || para.priorType == PriorType.LAPLACE_PRIOR
                || para.priorType == PriorType.TRUNC_NORMAL_PRIOR
                || (para.priorType == PriorType.GAMMA_PRIOR && para.gammaAlpha > 1)
                || (para.priorType == PriorType.GAMMA_PRIOR && para.gammaOffset > 0)
                || (para.priorType == PriorType.UNIFORM_PRIOR && para.uniformLower > 0 && para.uniformUpper < Double.POSITIVE_INFINITY)
                || (para.priorType == PriorType.EXPONENTIAL_PRIOR && para.exponentialOffset > 0)) {

            return true;
        } else {
            return false;
        }
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
        this.rateOptionClockModel = FixRateType.RElATIVE_TO;
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

    //+++++++++++++++++++++++ Validation ++++++++++++++++++++++++++++++++
    // true => valid, false => warning message 
    public boolean validateFixMeanRate(JCheckBox fixedMeanRateCheck) {
        return !(fixedMeanRateCheck.isSelected() && options.getPartitionClockModels().size() < 2);
    }
    
    public boolean validateRelativeTo() {
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            if (!model.isEstimatedRate()) { // fixed
                return true;
            }
        }        
        return false;
    }

}
