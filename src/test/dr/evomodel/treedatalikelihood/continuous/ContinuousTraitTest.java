/*
 * ContinuousTraitTest.java
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
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.*;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.io.StringReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class ContinuousTraitTest extends TraceCorrelationAssert {

    protected int dimTrait;
    protected int nTips;


    // Traits
    protected CompoundParameter traitParameter;

    // Diffusion
    protected MultivariateDiffusionModel diffusionModel;

    // Data Models
    protected ContinuousTraitPartialsProvider dataModel;
    ContinuousTraitPartialsProvider dataModelIntegrated;

    // Root
    protected ConjugateRootTraitPrior rootPrior;
    ConjugateRootTraitPrior rootPriorInf;
    ConjugateRootTraitPrior rootPriorIntegrated;

    // Factor
    protected MultivariateDiffusionModel diffusionModelFactor;
    protected IntegratedFactorAnalysisLikelihood dataModelFactor;
    protected Parameter rootMeanFactor;
    protected ConjugateRootTraitPrior rootPriorFactor;
    protected CompoundSymmetricMatrix precisionMatrixFactor;

    protected ContinuousRateTransformation rateTransformation;
    protected BranchRateModel rateModel;

    protected NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public ContinuousTraitTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        dimTrait = 3;

        format.setMaximumFractionDigits(5);

        // Tree
        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
        treeModel = createPrimateTreeModel();

        // Rates
        rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, false, false);
        rateModel = new DefaultBranchRateModel();

        // Data
        nTips = 6;
        Parameter[] dataTraits = new Parameter[6];
        dataTraits[0] = new Parameter.Default("human", new double[]{-1.0, 2.0, 3.0});
        dataTraits[1] = new Parameter.Default("chimp", new double[]{10.0, 12.0, 14.0});
        dataTraits[2] = new Parameter.Default("bonobo", new double[]{0.5, -2.0, 5.5});
        dataTraits[3] = new Parameter.Default("gorilla", new double[]{2.0, 5.0, -8.0});
        dataTraits[4] = new Parameter.Default("orangutan", new double[]{11.0, 1.0, -1.5});
        dataTraits[5] = new Parameter.Default("siamang", new double[]{1.0, 2.5, 4.0});
        traitParameter = new CompoundParameter("trait", dataTraits);

        List<Integer> missingIndices = new ArrayList<Integer>();
        traitParameter.setParameterValue(2, 0);
        missingIndices.add(3);
        missingIndices.add(4);
        missingIndices.add(5);
        missingIndices.add(7);

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
        Parameter rootMean = new Parameter.Default(new double[]{-1.0, -3.0, 2.5});
        Parameter rootSampleSize = new Parameter.Default(10.0);
        rootPrior = new ConjugateRootTraitPrior(rootMean, rootSampleSize);

        // Root prior Inf
        Parameter rootMeanInf = new Parameter.Default(new double[]{-1.0, -3.0, 2.0});
        Parameter rootSampleSizeInf = new Parameter.Default(Double.POSITIVE_INFINITY);
        rootPriorInf = new ConjugateRootTraitPrior(rootMeanInf, rootSampleSizeInf);

        // Data Model
        dataModel = new ContinuousTraitDataModel("dataModel",
                traitParameter,
                missingIndices, true,
                dimTrait, precisionType);

        //// Factor Model //// *****************************************************************************************
        // Diffusion
        Parameter offDiagonal = new Parameter.Default(new double[]{0.08164966});
        offDiagonal.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        Parameter diagonal = new Parameter.Default(new double[]{1.0, 1.5});

        precisionMatrixFactor = new CompoundSymmetricMatrix(diagonal, offDiagonal, true, true);
        diffusionModelFactor = new MultivariateDiffusionModel(precisionMatrixFactor);

        // Root prior
        rootMeanFactor = new Parameter.Default(new double[]{-1.0, 2.0});
        Parameter rootSampleSizeFac = new Parameter.Default(10.0);
        rootPriorFactor = new ConjugateRootTraitPrior(rootMeanFactor, rootSampleSizeFac);

        // Error model
        Parameter factorPrecisionParameters = new Parameter.Default("factorPrecision", new double[]{1.0, 5.0, 0.5});

        // Loadings
        Parameter[] loadingsParameters = new Parameter[2];
        loadingsParameters[0] = new Parameter.Default(new double[]{1.0, 2.0, 3.0});
        loadingsParameters[1] = new Parameter.Default(new double[]{0.0, 0.5, 1.0});
        MatrixParameterInterface loadingsMatrixParameters = new MatrixParameter("loadings", loadingsParameters);

        dataModelFactor = new IntegratedFactorAnalysisLikelihood("dataModelFactors",
                traitParameter,
                missingIndices,
                loadingsMatrixParameters,
                factorPrecisionParameters, 0.0, null,
                IntegratedFactorAnalysisLikelihood.CacheProvider.NO_CACHE);

        //// Integrated Process //// ***********************************************************************************
        // Data Model
        dataModelIntegrated = new IntegratedProcessTraitDataModel("dataModelIntegrated",
                traitParameter,
                missingIndices, true,
                dimTrait, precisionType);

        // Root prior
        String sI = "<beast>\n" +
                "    <conjugateRootPrior>\n" +
                "        <meanParameter>\n" +
                "            <parameter id=\"meanRoot\"  value=\"0.0 1.2 -0.5 -1.0 -3.0 2.5\"/>\n" +
                "        </meanParameter>\n" +
                "        <priorSampleSize>\n" +
                "            <parameter id=\"sampleSizeRoot\" value=\"10.0\"/>\n" +
                "        </priorSampleSize>\n" +
                "    </conjugateRootPrior>\n" +
                "</beast>";
        XMLParser parserI = new XMLParser(true, true, true, null);
        parserI.addXMLObjectParser(new AttributeParser());
        parserI.addXMLObjectParser(new ParameterParser());
        parserI.parse(new StringReader(sI), true);
        rootPriorIntegrated = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(parserI.getRoot(), 2 * dimTrait);
    }

    public void testDummy() {
    }
}


