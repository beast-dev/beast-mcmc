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

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.Vector;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.io.StringReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import static test.dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegateTest.getConditionalSimulations;
import static test.dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegateTest.getLogDatumLikelihood;

/**
 * @author Paul Bastide
 */

public class RepeatedMeasureFactorTest extends TraceCorrelationAssert {

    private int dimTrait;
    private int nTips;

    // Diffusion
    private MultivariateDiffusionModel diffusionModel;
    private ConjugateRootTraitPrior rootPrior;

    // Data Models
    private RepeatedMeasuresTraitDataModel dataModelRepeatedMeasures;
    private IntegratedFactorAnalysisLikelihood dataModelFactor;

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public RepeatedMeasureFactorTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        dimTrait = 3;

        format.setMaximumFractionDigits(5);

        // Tree
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();

        // Data
        nTips = 6;
        Parameter[] dataTraits = new Parameter[6];
        dataTraits[0] = new Parameter.Default("human", new double[]{-1.0, 2.0, 3.0});
        dataTraits[1] = new Parameter.Default("chimp", new double[]{10.0, 12.0, 14.0});
        dataTraits[2] = new Parameter.Default("bonobo", new double[]{0.5, -2.0, 5.5});
        dataTraits[3] = new Parameter.Default("gorilla", new double[]{2.0, 5.0, -8.0});
        dataTraits[4] = new Parameter.Default("orangutan", new double[]{11.0, 1.0, -1.5});
        dataTraits[5] = new Parameter.Default("siamang", new double[]{1.0, 2.5, 4.0});
        CompoundParameter traitParameter = new CompoundParameter("trait", dataTraits);

        List<Integer> missingIndices = new ArrayList<Integer>();
        traitParameter.setParameterValue(2, 0);
//        missingIndices.add(3);
//        missingIndices.add(4);
//        missingIndices.add(5);
//        missingIndices.add(7);

        //// Standard Model //// ***************************************************************************************

        // Diffusion
        Parameter[] precisionParameters = new Parameter[dimTrait];
        precisionParameters[0] = new Parameter.Default(new double[]{1.0, 0.1, 0.2});
        precisionParameters[1] = new Parameter.Default(new double[]{0.1, 2.0, 0.0});
        precisionParameters[2] = new Parameter.Default(new double[]{0.2, 0.0, 3.0});
        MatrixParameterInterface diffusionPrecisionMatrixParameter
                = new MatrixParameter("precisionMatrix", precisionParameters);
        diffusionModel = new MultivariateDiffusionModel(diffusionPrecisionMatrixParameter);

        PrecisionType precisionType = PrecisionType.FULL;

        // Root prior
        String s = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"-1.0 -3.0 2.5\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"10.0\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parser = new XMLParser(true, true, true, null);
        parser.addXMLObjectParser(new AttributeParser());
        parser.addXMLObjectParser(new ParameterParser());
        parser.parse(new StringReader(s), true);
        rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parser.getRoot(), dimTrait);

        // Error model
        Parameter[] samplingPrecision = new Parameter[dimTrait];
        samplingPrecision[0] = new Parameter.Default(new double[]{0.1, 0.0, 0.0});
        samplingPrecision[1] = new Parameter.Default(new double[]{0.0, 0.2, 0.0});
        samplingPrecision[2] = new Parameter.Default(new double[]{0.0, 0.0, 0.3});
        MatrixParameterInterface samplingPrecisionParameter
                = new MatrixParameter("samplingPrecisionMatrix", samplingPrecision);
        Parameter samplingPrecisionDiagonal = new Parameter.Default("samplingPrecisionMatrix", new double[]{0.1, 0.2, 0.3});

        //// Factor Model //// *****************************************************************************************
        // Loadings
        Parameter[] loadingsParameters = new Parameter[3];
        loadingsParameters[0] = new Parameter.Default(new double[]{1.0, 0.0, 0.0});
        loadingsParameters[1] = new Parameter.Default(new double[]{0.0, 1.0, 0.0});
        loadingsParameters[2] = new Parameter.Default(new double[]{0.0, 0.0, 1.0});
        MatrixParameterInterface loadingsMatrixParameters = new MatrixParameter("loadings", loadingsParameters);

        dataModelFactor = new IntegratedFactorAnalysisLikelihood("dataModelFactors",
                traitParameter,
                missingIndices,
                loadingsMatrixParameters,
                samplingPrecisionDiagonal, 0.0, null);


        //// Repeated Measures Model //// ******************************************************************************
        dataModelRepeatedMeasures = new RepeatedMeasuresTraitDataModel("dataModelRepeatedMeasures",
                traitParameter,
                missingIndices,
//                new boolean[3],
                true,
                dimTrait,
                samplingPrecisionParameter);

    }

    public void testLikelihoodBM() {
        System.out.println("\nTest Likelihood using vanilla BM:");

        // Diffusion
        DiffusionProcessDelegate diffusionProcessDelegate
                = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

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

        System.out.println("likelihoodBMRepMea: " + format.format(logDatumLikelihoodFactor));

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
    }

    public void testLikelihoodOU() {
        System.out.println("\nTest Likelihood using full OU:");

        // Diffusion
        List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        Parameter[] strengthOfSelectionParameters = new Parameter[3];
        strengthOfSelectionParameters[0] = new Parameter.Default(new double[]{0.5, 0.2, 0.0});
        strengthOfSelectionParameters[1] = new Parameter.Default(new double[]{0.2, 1.0, 0.1});
        strengthOfSelectionParameters[2] = new Parameter.Default(new double[]{0.0, 0.1, 10.5});
        MatrixParameter strengthOfSelectionMatrixParam
                = new MatrixParameter("strengthOfSelectionMatrix", strengthOfSelectionParameters);

        DiffusionProcessDelegate diffusionProcessDelegate
                = new OUDiffusionModelDelegate(treeModel, diffusionModel,
                optimalTraitsModels, new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        // Rates
        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        BranchRateModel rateModel = new DefaultBranchRateModel();

        //// Factor Model //// *****************************************************************************************
        // CDL
        ContinuousDataLikelihoodDelegate likelihoodDelegateFactors = new ContinuousDataLikelihoodDelegate(treeModel,
                diffusionProcessDelegate, dataModelFactor, rootPrior,
                rateTransformation, rateModel, true);

        dataModelFactor.setLikelihoodDelegate(likelihoodDelegateFactors);

        // Likelihood Computation
        TreeDataLikelihood dataLikelihoodFactors
                = new TreeDataLikelihood(likelihoodDelegateFactors, treeModel, rateModel);

        String sf = dataModelFactor.getReport();
        int indLikBegF = sf.indexOf("logMultiVariateNormalDensity = ") + 31;
        int indLikEndF = sf.indexOf("\n", indLikBegF);
        char[] logDatumLikelihoodCharF = new char[indLikEndF - indLikBegF + 1];
        sf.getChars(indLikBegF, indLikEndF, logDatumLikelihoodCharF, 0);
        double logDatumLikelihoodFactor = Double.parseDouble(String.valueOf(logDatumLikelihoodCharF));

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

        String s = dataLikelihoodRepMea.getReport();
        int indLikBeg = s.indexOf("logDatumLikelihood:") + 20;
        int indLikEnd = s.indexOf("\n", indLikBeg);
        char[] logDatumLikelihoodChar = new char[indLikEnd - indLikBeg + 1];
        s.getChars(indLikBeg, indLikEnd, logDatumLikelihoodChar, 0);
        double logDatumLikelihoodRepMea = Double.parseDouble(String.valueOf(logDatumLikelihoodChar));

        double likelihoodRepMeaDiffusion = dataLikelihoodRepMea.getLogLikelihood();

        assertEquals("likelihoodOURepMea",
                format.format(logDatumLikelihoodRepMea),
                format.format(likelihoodRepMeaDiffusion));

        System.out.println("likelihoodOURepMea: " + format.format(logDatumLikelihoodFactor));

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
    }
}
