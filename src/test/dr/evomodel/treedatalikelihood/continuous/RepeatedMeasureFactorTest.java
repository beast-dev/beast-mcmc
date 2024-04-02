/*
 * RepeatedMeasureFactorTest.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;

import java.util.ArrayList;
import java.util.List;

import static test.dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegateTest.getConditionalSimulations;
import static test.dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegateTest.getLogDatumLikelihood;

/**
 * @author Paul Bastide
 */

public class RepeatedMeasureFactorTest extends ContinuousTraitTest {

    // Data Models
    private RepeatedMeasuresTraitDataModel dataModelRepeatedMeasures;
    private RepeatedMeasuresTraitDataModel dataModelRepeatedMeasuresFull;

    public RepeatedMeasureFactorTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        dimTrait = 6;

        // Diffusion Model
        Parameter offDiagonal = new Parameter.Default(new double[]{0.12, -0.13, 0.14, -0.15, 0.16,
                -0.12, 0.13, -0.14, 0.15,
                0.12, -0.13, 0.14,
                -0.12, 0.13,
                0.12});
        offDiagonal.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, dimTrait * (dimTrait - 1) / 2));
        Parameter diagonal = new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6});
        diffusionModel = new MultivariateDiffusionModel(new CompoundSymmetricMatrix(diagonal, offDiagonal, true, true));

        // Root prior
        Parameter rootMean = new Parameter.Default(new double[]{-1.0, -3.0, 2.5, 1.0, -1.0, 0.0});
        Parameter rootSampleSize = new Parameter.Default(10.0);
        rootPrior = new ConjugateRootTraitPrior(rootMean, rootSampleSize);

        // Data
        Parameter[] dataTraits = new Parameter[6];
        dataTraits[0] = new Parameter.Default("human", new double[]{-1.0, 2.0, 3.0, 4.0, 5.0, -6.0});
        dataTraits[1] = new Parameter.Default("chimp", new double[]{10.0, 12.0, 14.0, 16.0, 18.0, 20.0});
        dataTraits[2] = new Parameter.Default("bonobo", new double[]{0.5, -2.0, 5.5, -5.2, 3.1, 1.1});
        dataTraits[3] = new Parameter.Default("gorilla", new double[]{2.0, 5.0, -8.0, -4.0, 3.2, 3.4});
        dataTraits[4] = new Parameter.Default("orangutan", new double[]{11.0, 1.0, -1.5, 2.4, -4.2, 6.0});
        dataTraits[5] = new Parameter.Default("siamang", new double[]{1.0, 2.5, 4.0, 4.0, -5.2, 1.0});
        CompoundParameter traitParameter = new CompoundParameter("trait", dataTraits);

        boolean[] missingIndicators = new boolean[traitParameter.getDimension()];
        traitParameter.setParameterValue(2, 0);
        missingIndicators[6] = true;
        missingIndicators[7] = true;
        missingIndicators[8] = true;
        missingIndicators[9] = true;
        missingIndicators[10] = true;
        missingIndicators[11] = true;
        missingIndicators[13] = true;
        missingIndicators[15] = true;
        missingIndicators[25] = true;
        missingIndicators[29] = true;

        // Error model Diagonal
        Parameter[] samplingPrecision = new Parameter[dimTrait];
        samplingPrecision[0] = new Parameter.Default(new double[]{0.1, 0.0, 0.0, 0.0, 0.0, 0.0});
        samplingPrecision[1] = new Parameter.Default(new double[]{0.0, 0.2, 0.0, 0.0, 0.0, 0.0});
        samplingPrecision[2] = new Parameter.Default(new double[]{0.0, 0.0, 0.3, 0.0, 0.0, 0.0});
        samplingPrecision[3] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.4, 0.0, 0.0});
        samplingPrecision[4] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.0, 0.5, 0.0});
        samplingPrecision[5] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.6});
        MatrixParameterInterface samplingPrecisionParameter
                = new MatrixParameter("samplingPrecisionMatrix", samplingPrecision);
        Parameter samplingPrecisionDiagonal = new Parameter.Default("samplingPrecisionMatrix", new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6});

        // Error model Full
        Parameter[] samplingPrecisionFull = new Parameter[dimTrait];
        samplingPrecision[0] = new Parameter.Default(new double[]{0.5, 0.0, 0.1, 0.1, 0.0, 0.0});
        samplingPrecision[1] = new Parameter.Default(new double[]{0.0, 0.2, 0.0, 0.0, 0.0, 0.0});
        samplingPrecision[2] = new Parameter.Default(new double[]{0.1, 0.0, 0.3, 0.0, 0.0, 0.2});
        samplingPrecision[3] = new Parameter.Default(new double[]{0.1, 0.0, 0.0, 0.4, 0.15, 0.0});
        samplingPrecision[4] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.15, 0.5, 0.0});
        samplingPrecision[5] = new Parameter.Default(new double[]{0.0, 0.0, 0.2, 0.0, 0.0, 0.6});
        MatrixParameterInterface samplingPrecisionParameterFull
                = new MatrixParameter("samplingPrecisionMatrix", samplingPrecision);

        //// Factor Model //// *****************************************************************************************
        // Loadings
        Parameter[] loadingsParameters = new Parameter[dimTrait];
        loadingsParameters[0] = new Parameter.Default(new double[]{1.0, 0.0, 0.0, 0.0, 0.0, 0.0});
        loadingsParameters[1] = new Parameter.Default(new double[]{0.0, 1.0, 0.0, 0.0, 0.0, 0.0});
        loadingsParameters[2] = new Parameter.Default(new double[]{0.0, 0.0, 1.0, 0.0, 0.0, 0.0});
        loadingsParameters[3] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 1.0, 0.0, 0.0});
        loadingsParameters[4] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.0, 1.0, 0.0});
        loadingsParameters[5] = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
        MatrixParameterInterface loadingsMatrixParameters = new MatrixParameter("loadings", loadingsParameters);

        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndicators,
                true,
                6,
                PrecisionType.FULL
        );

        dataModelFactor = new IntegratedFactorAnalysisLikelihood("dataModelFactors",
                traitParameter,
                missingIndicators,
                loadingsMatrixParameters,
                samplingPrecisionDiagonal, 0.0, null,
                IntegratedFactorAnalysisLikelihood.CacheProvider.NO_CACHE);


        //// Repeated Measures Model //// ******************************************************************************
        dataModelRepeatedMeasures = new RepeatedMeasuresTraitDataModel("dataModelRepeatedMeasures",
                dataModel,
                traitParameter,
                missingIndicators,
//                new boolean[3],
                true,
                dimTrait,
                1,
                samplingPrecisionParameter,
                PrecisionType.FULL);

        dataModelRepeatedMeasuresFull = new RepeatedMeasuresTraitDataModel("dataModelRepeatedMeasures",
                dataModel,
                traitParameter,
                missingIndicators,
                true,
                dimTrait,
                1,
                samplingPrecisionParameterFull,
                PrecisionType.FULL);

    }

    public void testLikelihoodBM() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        //// Factor Model //// *****************************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPrior,
                rateTransformation, rateModel, true);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        double logDatumLikelihoodFactor = getLogDatumLikelihood(dataModelFactor);

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();

        assertEquals("likelihoodBMFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        System.out.println("likelihoodBMFactor: " + format.format(logDatumLikelihoodFactor));

        // Simulation
        MathUtils.setSeed(17890826);
        double[] traitsFactors = getConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModel, dataModelFactor, rootPrior, treeModel, rateTransformation);
        System.err.println(new Vector(traitsFactors));

        //// Repeated Measures //// ************************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateRepMea = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior,
                rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodRepMea
                = new TreeDataLikelihood(likelihoodDelegateRepMea, treeModel, rateModel);

        double logDatumLikelihoodRepMea = getLogDatumLikelihood(dataLikelihoodRepMea);

        double likelihoodRepMeaDiffusion = dataLikelihoodRepMea.getLogLikelihood();

        assertEquals("likelihoodBMRepMea",
                format.format(logDatumLikelihoodRepMea),
                format.format(likelihoodRepMeaDiffusion));

        System.out.println("likelihoodBMRepMea: " + format.format(logDatumLikelihoodRepMea));

        // Simulation
        MathUtils.setSeed(17890826);
        double[] traitsRepMea = getConditionalSimulations(dataLikelihoodRepMea, likelihoodDelegateRepMea, diffusionModel, dataModelRepeatedMeasures, rootPrior, treeModel, rateTransformation);
        System.err.println(new Vector(traitsRepMea));

        //// Equal ? //// **********************************************************************************************
        assertEquals("likelihoodBMRepFactor",
                format.format(likelihoodFactorData + likelihoodFactorDiffusion),
                format.format(likelihoodRepMeaDiffusion));

        for (int i = 0; i < traitsFactors.length; i++) {
            assertEquals(format.format(traitsRepMea[i]), format.format(traitsFactors[i]));
        }

        //// Repeated Measures Full //// *******************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateRepMeaFull = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelRepeatedMeasuresFull, rootPrior,
                rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodRepMeaFull
                = new TreeDataLikelihood(likelihoodDelegateRepMeaFull, treeModel, rateModel);

        double logDatumLikelihoodRepMeaFull = getLogDatumLikelihood(dataLikelihoodRepMeaFull);

        double likelihoodRepMeaDiffusionFull = dataLikelihoodRepMeaFull.getLogLikelihood();

        assertEquals("likelihoodBMRepMea",
                format.format(logDatumLikelihoodRepMeaFull),
                format.format(likelihoodRepMeaDiffusionFull));

        System.out.println("likelihoodBMRepMeaFull: " + format.format(logDatumLikelihoodRepMeaFull));
    }

    public void testLikelihoodOU() {
        System.out.println("\nTest Likelihood using full OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.4", new double[]{10.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.5", new double[]{20.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.6", new double[]{-20.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[6];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{1.0, 0.1, 0.0, 0.0, 0.5, 2.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.1, 10., 0.0, 0.0, 0.0, 0.0});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.0, 20., 0.3, 0.0, 0.0});
        strengthOfSelectionParameters[3] = new Parameter.Default(new double[]{0.0, 0.0, 0.3, 30., 3.0, 0.0});
        strengthOfSelectionParameters[4] = new Parameter.Default(new double[]{1.0, 0.0, 0.0, 3.0, 40., 0.0});
        strengthOfSelectionParameters[5] = new Parameter.Default(new double[]{0.0, 0.0, 0.5, 0.0, 0.0, 50.});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        //// Factor Model //// *****************************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPrior,
                rateTransformation, rateModel, true);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        double logDatumLikelihoodFactor = getLogDatumLikelihood(dataModelFactor);

        double likelihoodFactorData = dataLikelihoodFactors.getLogLikelihood();
        double likelihoodFactorDiffusion = dataModelFactor.getLogLikelihood();

        assertEquals("likelihoodOUFactor",
                format.format(logDatumLikelihoodFactor),
                format.format(likelihoodFactorData + likelihoodFactorDiffusion));

        System.out.println("likelihoodOUFactor: " + format.format(logDatumLikelihoodFactor));

        // Simulation
        MathUtils.setSeed(17890826);
        double[] traitsFactors = getConditionalSimulations(dataLikelihoodFactors, likelihoodDelegateFactors, diffusionModel, dataModelFactor, rootPrior, treeModel, rateTransformation);
        System.err.println(new Vector(traitsFactors));

        //// Repeated Measures //// ************************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateRepMea = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelRepeatedMeasures, rootPrior,
                rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodRepMea
                = new TreeDataLikelihood(likelihoodDelegateRepMea, treeModel, rateModel);

        double logDatumLikelihoodRepMea = getLogDatumLikelihood(dataLikelihoodRepMea);

        double likelihoodRepMeaDiffusion = dataLikelihoodRepMea.getLogLikelihood();

        assertEquals("likelihoodOURepMea",
                format.format(logDatumLikelihoodRepMea),
                format.format(likelihoodRepMeaDiffusion));

        System.out.println("likelihoodOURepMea: " + format.format(logDatumLikelihoodRepMea));

        // Simulation
        MathUtils.setSeed(17890826);
        double[] traitsRepMea = getConditionalSimulations(dataLikelihoodRepMea, likelihoodDelegateRepMea, diffusionModel, dataModelRepeatedMeasures, rootPrior, treeModel, rateTransformation);
        System.err.println(new Vector(traitsRepMea));

        //// Equal ? //// **********************************************************************************************
        assertEquals("likelihoodOURepFactor",
                format.format(likelihoodFactorData + likelihoodFactorDiffusion),
                format.format(likelihoodRepMeaDiffusion));

        for (int i = 0; i < traitsFactors.length; i++) {
            assertEquals(format.format(traitsRepMea[i]), format.format(traitsFactors[i]));
        }

        //// Repeated Measures Full //// *******************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateRepMeaFull = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelRepeatedMeasuresFull, rootPrior,
                rateTransformation, rateModel, true);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodRepMeaFull
                = new TreeDataLikelihood(likelihoodDelegateRepMeaFull, treeModel, rateModel);

        double logDatumLikelihoodRepMeaFull = getLogDatumLikelihood(dataLikelihoodRepMeaFull);

        double likelihoodRepMeaDiffusionFull = dataLikelihoodRepMeaFull.getLogLikelihood();

        assertEquals("likelihoodBMRepMea",
                format.format(logDatumLikelihoodRepMeaFull),
                format.format(likelihoodRepMeaDiffusionFull));

        System.out.println("likelihoodBMRepMeaFull: " + format.format(logDatumLikelihoodRepMeaFull));
    }
}
