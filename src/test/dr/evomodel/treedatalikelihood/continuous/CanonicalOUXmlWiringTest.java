/*
 * CanonicalOUXmlWiringTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.adapter.CanonicalOUMessagePasserComputer;
import dr.evomodel.treedatalikelihood.continuous.adapter.CanonicalTipObservationAdapter;
import dr.evomodel.treedatalikelihood.hmc.CanonicalMeanParameterGradient;
import dr.evomodel.treedatalikelihood.hmc.CanonicalPrecisionGradient;
import dr.evomodel.treedatalikelihood.hmc.CanonicalSelectionParameterGradient;
import dr.inference.model.CachedMatrixInverse;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.CompoundParameter;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

public class CanonicalOUXmlWiringTest extends ContinuousTraitTest {

    private static final double FD_DELTA = 1.0e-5;
    private static final double TOL = 1.0e-4;

    private static final double ORTHOGONAL_SCALAR = 1.3;
    private static final double ORTHOGONAL_RHO = 0.85;
    private static final double ORTHOGONAL_THETA = 0.25;
    private static final double ORTHOGONAL_T = -0.08;
    private static final double[] ORTHOGONAL_ANGLES = {0.2, -0.15, 0.1};
    private static final double[] OPTIMA = {1.0, -0.5, 2.0};

    public CanonicalOUXmlWiringTest(final String name) {
        super(name);
    }

    public void testCanonicalDelegateLikelihoodMatchesStandaloneComputerWithMissingData() {
        final XmlStyleCanonicalOUSetup setup = buildCanonicalSetup("xmlCanonLike");
        final double delegateLogLikelihood = setup.treeDataLikelihood.getLogLikelihood();
        final double standaloneLogLikelihood = setup.standaloneComputer.computeLogLikelihood();
        assertEquals("canonical XML-style delegate log-likelihood", standaloneLogLikelihood, delegateLogLikelihood, 1.0e-8);
    }

    public void testCanonicalDelegateOrthogonalBlockSelectionGradientMatchesFiniteDifferenceWithMissingData() {
        final XmlStyleCanonicalOUSetup setup = buildCanonicalSetup("xmlCanonGrad");
        final CanonicalSelectionParameterGradient gradientProvider =
                new CanonicalSelectionParameterGradient(
                        setup.treeDataLikelihood,
                        setup.delegate,
                        setup.nativeParameter,
                        setup.blockSelection);

        final double[] analytic = gradientProvider.getGradientLogDensity();
        final double[] numeric = new double[setup.nativeParameter.getDimension()];

        for (int i = 0; i < numeric.length; i++) {
            final double original = setup.nativeParameter.getParameterValue(i);

            setup.nativeParameter.setParameterValue(i, original + FD_DELTA);
            setup.treeDataLikelihood.makeDirty();
            final double plus = setup.treeDataLikelihood.getLogLikelihood();

            setup.nativeParameter.setParameterValue(i, original - FD_DELTA);
            setup.treeDataLikelihood.makeDirty();
            final double minus = setup.treeDataLikelihood.getLogLikelihood();

            setup.nativeParameter.setParameterValue(i, original);
            numeric[i] = (plus - minus) / (2.0 * FD_DELTA);
        }
        setup.treeDataLikelihood.makeDirty();

        for (int i = 0; i < numeric.length; i++) {
            assertEquals("canonical XML-style native selection gradient[" + i + "]",
                    numeric[i], analytic[i], TOL);
        }
    }

    public void testCanonicalDelegateFixedRootMeanGradientMatchesFiniteDifference() {
        final XmlStyleCanonicalMeanSetup setup = buildCanonicalFixedRootMeanSetup("xmlCanonMeanFixed");
        final CanonicalMeanParameterGradient gradientProvider =
                new CanonicalMeanParameterGradient(
                        setup.treeDataLikelihood,
                        setup.delegate,
                        setup.meanParameter,
                        true,
                        true);

        final double[] analytic = gradientProvider.getGradientLogDensity();
        final double[] numeric = new double[setup.meanParameter.getDimension()];

        for (int i = 0; i < numeric.length; i++) {
            final double original = setup.meanParameter.getParameterValue(i);

            setup.meanParameter.setParameterValue(i, original + FD_DELTA);
            setup.treeDataLikelihood.makeDirty();
            final double plus = setup.treeDataLikelihood.getLogLikelihood();

            setup.meanParameter.setParameterValue(i, original - FD_DELTA);
            setup.treeDataLikelihood.makeDirty();
            final double minus = setup.treeDataLikelihood.getLogLikelihood();

            setup.meanParameter.setParameterValue(i, original);
            numeric[i] = (plus - minus) / (2.0 * FD_DELTA);
        }
        setup.treeDataLikelihood.makeDirty();

        for (int i = 0; i < numeric.length; i++) {
            assertEquals("canonical XML-style fixed-root mean gradient[" + i + "]",
                    numeric[i], analytic[i], TOL);
        }
    }

    public void testCanonicalDelegatePrecisionGradientMatchesFiniteDifference() {
        final XmlStyleCanonicalPrecisionSetup setup = buildCanonicalPrecisionSetup("xmlCanonPrecision");
        final CanonicalPrecisionGradient gradientProvider =
                new CanonicalPrecisionGradient(
                        setup.treeDataLikelihood,
                        setup.delegate,
                        setup.precisionMatrix);

        final double[] analytic = gradientProvider.getGradientLogDensity();
        final double[] numeric = new double[gradientProvider.getDimension()];
        final Parameter parameter = gradientProvider.getParameter();

        for (int i = 0; i < numeric.length; i++) {
            final double original = parameter.getParameterValue(i);

            parameter.setParameterValue(i, original + FD_DELTA);
            setup.treeDataLikelihood.makeDirty();
            final double plus = setup.treeDataLikelihood.getLogLikelihood();

            parameter.setParameterValue(i, original - FD_DELTA);
            setup.treeDataLikelihood.makeDirty();
            final double minus = setup.treeDataLikelihood.getLogLikelihood();

            parameter.setParameterValue(i, original);
            numeric[i] = (plus - minus) / (2.0 * FD_DELTA);
        }
        setup.treeDataLikelihood.makeDirty();

        for (int i = 0; i < numeric.length; i++) {
            assertEquals("canonical XML-style precision gradient[" + i + "]",
                    numeric[i], analytic[i], TOL);
        }
    }

    private XmlStyleCanonicalOUSetup buildCanonicalSetup(final String tag) {
        final Parameter angles = new Parameter.Default("angles." + tag, ORTHOGONAL_ANGLES.clone());
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation." + tag, dimTrait, angles);
        final Parameter scalar = new Parameter.Default("scalar." + tag, new double[]{ORTHOGONAL_SCALAR});
        final Parameter rho = new Parameter.Default("rho." + tag, new double[]{ORTHOGONAL_RHO});
        final Parameter theta = new Parameter.Default("theta." + tag, new double[]{ORTHOGONAL_THETA});
        final Parameter t = new Parameter.Default("t." + tag, new double[]{ORTHOGONAL_T});

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "selection." + tag,
                        rotation,
                        scalar,
                        rho,
                        theta,
                        t);

        final List<BranchRateModel> optima = new ArrayList<BranchRateModel>(dimTrait);
        for (int i = 0; i < dimTrait; i++) {
            optima.add(new StrictClockBranchRates(
                    new Parameter.Default("optimum." + tag + "." + i, new double[]{OPTIMA[i]})));
        }

        final OUDiffusionModelDelegate diffusionDelegate =
                new OUDiffusionModelDelegate(
                        treeModel,
                        diffusionModel,
                        optima,
                        new MultivariateElasticModel(blockSelection));

        final ContinuousDataLikelihoodDelegate delegate =
                new ContinuousDataLikelihoodDelegate(
                        treeModel,
                        diffusionDelegate,
                        dataModel,
                        rootPrior,
                        rateTransformation,
                        rateModel,
                        false,
                        false);
        delegate.enableCanonicalOULikelihood();

        final TreeDataLikelihood treeDataLikelihood =
                new TreeDataLikelihood(delegate, treeModel, rateModel);

        final CanonicalOUMessagePasserComputer standaloneComputer =
                new CanonicalOUMessagePasserComputer(
                        treeModel,
                        diffusionDelegate.getElasticModel(),
                        diffusionModel,
                        CanonicalTipObservationAdapter.extractTipObservations(treeModel, dataModel, dimTrait),
                        rootPrior,
                        diffusionDelegate.getCanonicalStationaryMeanParameter(),
                        rateModel,
                        rateTransformation);

        return new XmlStyleCanonicalOUSetup(
                blockSelection,
                blockSelection.getParameter(),
                delegate,
                treeDataLikelihood,
                standaloneComputer);
    }

    private XmlStyleCanonicalMeanSetup buildCanonicalFixedRootMeanSetup(final String tag) {
        final Parameter angles = new Parameter.Default("angles." + tag, ORTHOGONAL_ANGLES.clone());
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation." + tag, dimTrait, angles);
        final Parameter scalar = new Parameter.Default("scalar." + tag, new double[]{ORTHOGONAL_SCALAR});
        final Parameter rho = new Parameter.Default("rho." + tag, new double[]{ORTHOGONAL_RHO});
        final Parameter theta = new Parameter.Default("theta." + tag, new double[]{ORTHOGONAL_THETA});
        final Parameter t = new Parameter.Default("t." + tag, new double[]{ORTHOGONAL_T});

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "selection." + tag,
                        rotation,
                        scalar,
                        rho,
                        theta,
                        t);

        final List<BranchRateModel> optima = new ArrayList<BranchRateModel>(dimTrait);
        final CompoundParameter meanParameter = new CompoundParameter("mean." + tag);
        for (int i = 0; i < dimTrait; i++) {
            final Parameter optimum = new Parameter.Default("optimum." + tag + "." + i, new double[]{OPTIMA[i]});
            meanParameter.addParameter(optimum);
            optima.add(new StrictClockBranchRates(optimum));
        }

        final OUDiffusionModelDelegate diffusionDelegate =
                new OUDiffusionModelDelegate(
                        treeModel,
                        diffusionModel,
                        optima,
                        new MultivariateElasticModel(blockSelection));

        final ConjugateRootTraitPrior fixedRootPrior =
                new ConjugateRootTraitPrior(
                        meanParameter,
                        new Parameter.Default(Double.POSITIVE_INFINITY));

        final ContinuousDataLikelihoodDelegate delegate =
                new ContinuousDataLikelihoodDelegate(
                        treeModel,
                        diffusionDelegate,
                        dataModel,
                        fixedRootPrior,
                        rateTransformation,
                        rateModel,
                        false,
                        false);
        delegate.enableCanonicalOULikelihood();

        final TreeDataLikelihood treeDataLikelihood =
                new TreeDataLikelihood(delegate, treeModel, rateModel);

        return new XmlStyleCanonicalMeanSetup(meanParameter, delegate, treeDataLikelihood);
    }

    private XmlStyleCanonicalPrecisionSetup buildCanonicalPrecisionSetup(final String tag) {
        final Parameter angles = new Parameter.Default("angles." + tag, ORTHOGONAL_ANGLES.clone());
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation." + tag, dimTrait, angles);
        final Parameter scalar = new Parameter.Default("scalar." + tag, new double[]{ORTHOGONAL_SCALAR});
        final Parameter rho = new Parameter.Default("rho." + tag, new double[]{ORTHOGONAL_RHO});
        final Parameter theta = new Parameter.Default("theta." + tag, new double[]{ORTHOGONAL_THETA});
        final Parameter t = new Parameter.Default("t." + tag, new double[]{ORTHOGONAL_T});

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "selection." + tag,
                        rotation,
                        scalar,
                        rho,
                        theta,
                        t);

        final Parameter diagonal = new Parameter.Default("varianceDiagonal." + tag, new double[]{1.2, 0.9, 1.5});
        final Parameter offDiagonal = new Parameter.Default("varianceOffDiagonal." + tag, new double[]{0.15, -0.05, 0.08});
        final CompoundSymmetricMatrix varianceMatrix =
                new CompoundSymmetricMatrix(diagonal, offDiagonal, true, false);
        final CachedMatrixInverse precisionMatrix = new CachedMatrixInverse("precision." + tag, varianceMatrix);
        final dr.evomodel.continuous.MultivariateDiffusionModel localDiffusionModel =
                new dr.evomodel.continuous.MultivariateDiffusionModel(precisionMatrix);

        final List<BranchRateModel> optima = new ArrayList<BranchRateModel>(dimTrait);
        for (int i = 0; i < dimTrait; i++) {
            optima.add(new StrictClockBranchRates(
                    new Parameter.Default("optimum." + tag + "." + i, new double[]{OPTIMA[i]})));
        }

        final OUDiffusionModelDelegate diffusionDelegate =
                new OUDiffusionModelDelegate(
                        treeModel,
                        localDiffusionModel,
                        optima,
                        new MultivariateElasticModel(blockSelection));

        final ContinuousDataLikelihoodDelegate delegate =
                new ContinuousDataLikelihoodDelegate(
                        treeModel,
                        diffusionDelegate,
                        dataModel,
                        rootPrior,
                        rateTransformation,
                        rateModel,
                        false,
                        false);
        delegate.enableCanonicalOULikelihood();

        final TreeDataLikelihood treeDataLikelihood =
                new TreeDataLikelihood(delegate, treeModel, rateModel);

        return new XmlStyleCanonicalPrecisionSetup(precisionMatrix, delegate, treeDataLikelihood);
    }

    private static final class XmlStyleCanonicalOUSetup {
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection;
        final Parameter nativeParameter;
        final ContinuousDataLikelihoodDelegate delegate;
        final TreeDataLikelihood treeDataLikelihood;
        final CanonicalOUMessagePasserComputer standaloneComputer;

        private XmlStyleCanonicalOUSetup(final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection,
                                         final Parameter nativeParameter,
                                         final ContinuousDataLikelihoodDelegate delegate,
                                         final TreeDataLikelihood treeDataLikelihood,
                                         final CanonicalOUMessagePasserComputer standaloneComputer) {
            this.blockSelection = blockSelection;
            this.nativeParameter = nativeParameter;
            this.delegate = delegate;
            this.treeDataLikelihood = treeDataLikelihood;
            this.standaloneComputer = standaloneComputer;
        }
    }

    private static final class XmlStyleCanonicalMeanSetup {
        final CompoundParameter meanParameter;
        final ContinuousDataLikelihoodDelegate delegate;
        final TreeDataLikelihood treeDataLikelihood;

        private XmlStyleCanonicalMeanSetup(final CompoundParameter meanParameter,
                                           final ContinuousDataLikelihoodDelegate delegate,
                                           final TreeDataLikelihood treeDataLikelihood) {
            this.meanParameter = meanParameter;
            this.delegate = delegate;
            this.treeDataLikelihood = treeDataLikelihood;
        }
    }

    private static final class XmlStyleCanonicalPrecisionSetup {
        final MatrixParameterInterface precisionMatrix;
        final ContinuousDataLikelihoodDelegate delegate;
        final TreeDataLikelihood treeDataLikelihood;

        private XmlStyleCanonicalPrecisionSetup(final MatrixParameterInterface precisionMatrix,
                                                final ContinuousDataLikelihoodDelegate delegate,
                                                final TreeDataLikelihood treeDataLikelihood) {
            this.precisionMatrix = precisionMatrix;
            this.delegate = delegate;
            this.treeDataLikelihood = treeDataLikelihood;
        }
    }
}
