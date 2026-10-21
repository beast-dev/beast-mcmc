/*
 * ContinuousDataLikelihoodParser.java
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

package dr.evomodelxml.treedatalikelihood;

import dr.app.beauti.components.BeautiModelIDProvider;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUBranchProcessSpecification;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUOptions;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuousNew.ContinuousDataLikelihoodDelegateNew;
import dr.evomodel.treedatalikelihood.continuousNew.CtsCompatibility;
import dr.evomodel.treedatalikelihood.preorder.*;
import dr.evomodelxml.continuous.ContinuousTraitDataModelParser;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inferencexml.timeseries.OUSelectionChartParserHelper;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.getTipTraitName;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 */
public class ContinuousDataLikelihoodParser extends AbstractXMLObjectParser implements BeautiModelIDProvider {

    public static final String CONJUGATE_ROOT_PRIOR = AbstractMultivariateTraitLikelihood.CONJUGATE_ROOT_PRIOR;
    private static final String USE_TREE_LENGTH = AbstractMultivariateTraitLikelihood.USE_TREE_LENGTH;
    private static final String SCALE_BY_TIME = AbstractMultivariateTraitLikelihood.SCALE_BY_TIME;
    private static final String RECIPROCAL_RATES = AbstractMultivariateTraitLikelihood.RECIPROCAL_RATES;
    private static final String DRIFT_MODELS = AbstractMultivariateTraitLikelihood.DRIFT_MODELS;
    private static final String OPTIMAL_TRAITS = AbstractMultivariateTraitLikelihood.OPTIMAL_TRAITS;
    private static final String REPORT_MULTIVARIATE = AbstractMultivariateTraitLikelihood.REPORT_MULTIVARIATE;
    private static final String INTEGRATE_INTERNAL_TRAITS = AbstractMultivariateTraitLikelihood.INTEGRATE;

    private static final String RECONSTRUCT_TRAITS = "reconstructTraits";
    private static final String ALLOW_SINGULAR = "allowSingular";
    public static final String FORCE_FULL_PRECISION = "forceFullPrecision"; // TODO: maybe discourage using this
    private static final String FORCE_COMPLETELY_MISSING = "forceCompletelyMissing"; // TODO: maybe discourage using this
    private static final String FORCE_COMPLETELY_OBSERVED = "forceCompletelyObserved";

    private static final String FORCE_DRIFT = "forceDrift";
    private static final String FORCE_OU = "forceOU";
    private static final String IMPLEMENTATION = "implementation";
    private static final String IMPLEMENTATION_LEGACY = "legacy";
    private static final String IMPLEMENTATION_CANONICAL = "canonical";
    private static final String IMPLEMENTATION_CDI = "cdi";
    private static final String IMPLEMENTATION_CDI_MOMENT = "cdi-moment";
    private static final String IMPLEMENTATION_BM_CDI = "bm-cdi";
    private static final String IMPLEMENTATION_BM_CDI_MOMENT = "bm-cdi-moment";
    private static final String IMPLEMENTATION_BM_MOMENT = "bm-moment";
    private static final String IMPLEMENTATION_NATIVE_MOMENT_BM = "native-moment-bm";
    private static final String IMPLEMENTATION_CANONICAL_MESSAGE = "canonical-message";
    private static final String IMPLEMENTATION_OU_CANONICAL = "ou-canonical";
    private static final String IMPLEMENTATION_OU_CANONICAL_MESSAGE = "ou-canonical-message";
    private static final String IMPLEMENTATION_NATURAL_CANONICAL_BM = "natural-canonical-bm";
    private static final String IMPLEMENTATION_BM_NATURAL_CANONICAL = "bm-natural-canonical";

    private static final String STRENGTH_OF_SELECTION_MATRIX = "strengthOfSelectionMatrix";
    private static final String CANONICAL_OU_BRANCH_PROCESSES = "canonicalOUBranchProcesses";
    private static final String BRANCH_PROCESS_INDEX = "branchProcessIndex";
    private static final String PROCESS = "process";
    private static final String CANONICAL_OU_COORDINATES = CanonicalOUOptions.COORDINATES_ATTRIBUTE;

    public static final String CONTINUOUS_DATA_LIKELIHOOD = "traitDataLikelihood";

    public String getParserName() {
        return CONTINUOUS_DATA_LIKELIHOOD;
    }

    private enum DelegateProvider {
        OU,
        DRIFT,
        HOMOGENOUS
    }

    private enum ImplementationMode {
        LEGACY,
        CANONICAL_COMPATIBILITY,
        CDI_MOMENT_BM,
        MOMENT_BM,
        CANONICAL_MESSAGE_OU,
        NATURAL_CANONICAL_BM;

        static ImplementationMode parse(final String rawImplementation) throws XMLParseException {
            final String implementation = rawImplementation.trim().toLowerCase(Locale.ROOT);
            if (IMPLEMENTATION_LEGACY.equals(implementation)) {
                return LEGACY;
            }
            if (IMPLEMENTATION_CANONICAL.equals(implementation)) {
                return CANONICAL_COMPATIBILITY;
            }
            if (IMPLEMENTATION_CDI.equals(implementation)
                    || IMPLEMENTATION_CDI_MOMENT.equals(implementation)
                    || IMPLEMENTATION_BM_CDI.equals(implementation)
                    || IMPLEMENTATION_BM_CDI_MOMENT.equals(implementation)) {
                return CDI_MOMENT_BM;
            }
            if (IMPLEMENTATION_BM_MOMENT.equals(implementation)
                    || IMPLEMENTATION_NATIVE_MOMENT_BM.equals(implementation)) {
                return MOMENT_BM;
            }
            if (IMPLEMENTATION_CANONICAL_MESSAGE.equals(implementation)
                    || IMPLEMENTATION_OU_CANONICAL.equals(implementation)
                    || IMPLEMENTATION_OU_CANONICAL_MESSAGE.equals(implementation)) {
                return CANONICAL_MESSAGE_OU;
            }
            if (IMPLEMENTATION_NATURAL_CANONICAL_BM.equals(implementation)
                    || IMPLEMENTATION_BM_NATURAL_CANONICAL.equals(implementation)) {
                return NATURAL_CANONICAL_BM;
            }

            throw new XMLParseException("traitDataLikelihood implementation must be one of "
                    + quotedImplementationList() + ", not \"" + rawImplementation + "\".");
        }

        boolean requiresLikelihoodOnly() {
            return this != LEGACY;
        }
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // Begin Parse Tree and rates
        Tree treeModel = (Tree) xo.getChild(Tree.class);
        MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        boolean useTreeLength = xo.getAttribute(USE_TREE_LENGTH, false);
        boolean scaleByTime = xo.getAttribute(SCALE_BY_TIME, false);
        boolean reciprocalRates = xo.getAttribute(RECIPROCAL_RATES, false);
        final String implementation = xo.getAttribute(IMPLEMENTATION, IMPLEMENTATION_LEGACY);
        final ImplementationMode implementationMode = ImplementationMode.parse(implementation);
        final CanonicalOUOptions canonicalOUOptions = parseCanonicalOUOptions(xo);

        if (reciprocalRates) {
            throw new XMLParseException("Reciprocal rates are not yet implemented.");
        }

        if (rateModel == null) {
            rateModel = new DefaultBranchRateModel();
        }

        ContinuousRateTransformation rateTransformation = new ContinuousRateTransformation.Default(
                treeModel, scaleByTime, useTreeLength);
        // End Parse Tree and rates


        // Begin Parse Evolution Model
        List<BranchRateModel> driftModels = AbstractMultivariateTraitLikelihood.parseDriftModels(xo, diffusionModel);
        List<BranchRateModel> optimalTraitsModels = AbstractMultivariateTraitLikelihood.parseOptimalValuesModels(xo, diffusionModel);


        MultivariateElasticModel elasticModel = null;
        if (xo.hasChildNamed(STRENGTH_OF_SELECTION_MATRIX)) {
            XMLObject cxo = xo.getChild(STRENGTH_OF_SELECTION_MATRIX);
            MatrixParameterInterface strengthOfSelectionMatrixParam;
            strengthOfSelectionMatrixParam = (MatrixParameterInterface) cxo.getChild(MatrixParameterInterface.class);
            if (strengthOfSelectionMatrixParam != null) {
                elasticModel = new MultivariateElasticModel(strengthOfSelectionMatrixParam);
            }
        }
        final CanonicalOUBranchProcessSpecification branchProcessSpecification =
                parseCanonicalOUBranchProcesses(xo);
        if (branchProcessSpecification != null) {
            if (elasticModel != null) {
                throw new XMLParseException(CONTINUOUS_DATA_LIKELIHOOD + " cannot combine "
                        + "<" + STRENGTH_OF_SELECTION_MATRIX + "> with <"
                        + CANONICAL_OU_BRANCH_PROCESSES + ">.");
            }
            validateBranchProcessSelectionCharts(xo, branchProcessSpecification);
            elasticModel = branchProcessSpecification.getElasticModel(0);
        } else {
            OUSelectionChartParserHelper.validateSelectionChart(xo,
                    elasticModel == null ? null : elasticModel.getStrengthOfSelectionMatrixParameter(),
                    CONTINUOUS_DATA_LIKELIHOOD);
        }


        DiffusionProcessDelegate diffusionProcessDelegate;

        boolean forceOU = xo.getAttribute(FORCE_OU, false);
        boolean forceDrift = xo.getAttribute(FORCE_DRIFT, false);

        final DelegateProvider delegateProvider;

        if ((optimalTraitsModels != null && elasticModel != null) || forceOU) {
            delegateProvider = DelegateProvider.OU;
        } else {
            if (driftModels != null || forceDrift) {
                delegateProvider = DelegateProvider.DRIFT;
            } else {
                delegateProvider = DelegateProvider.HOMOGENOUS;
            }
        }


        final String traitName;

        ContinuousTraitPartialsProvider dataModel;

        if (xo.hasChildNamed(TreeTraitParserUtilities.TRAIT_PARAMETER)) {
            PrecisionType precisionType = null;
            if (delegateProvider != DelegateProvider.HOMOGENOUS) {
                precisionType = PrecisionType.FULL;
            }
            dataModel = ContinuousTraitDataModelParser.parseContinuousTraitDataModel(xo, precisionType);
        } else {  // Has ContinuousTraitPartialsProvider
            dataModel = (ContinuousTraitPartialsProvider) xo.getChild(ContinuousTraitPartialsProvider.class);

            if (delegateProvider != DelegateProvider.HOMOGENOUS && dataModel.getPrecisionType() != PrecisionType.FULL) {
                throw new RuntimeException("The data likelihood delegate is not 'HomogeneousDiffusionModelDelegate', " +
                        "but the data model with name '" + dataModel.getModelName() + "' has precision type '" +
                        dataModel.getPrecisionType().getTag() + "'.");
            }

            if (xo.hasChildNamed(TreeTraitParserUtilities.JITTER)) {
                System.err.println("Jitter is specified in " + xo.getAttribute(XMLObject.ID) + " but will be ignored, as a data model is provided. If a jitter is needed, please add it to the data model " + dataModel.getClass().toString() + ".");
            }
            if (xo.getAttribute(FORCE_FULL_PRECISION, false)) {
                System.err.println("Option " + FORCE_FULL_PRECISION + "is set to 'true' in " + xo.getAttribute(XMLObject.ID) + " but will be ignored, as a data model is provided. Please check the data model for this attribute.");
            }
            if (xo.getAttribute(FORCE_COMPLETELY_MISSING, false)) {
                System.err.println("Option " + FORCE_COMPLETELY_MISSING + "is set to 'true' in " + xo.getAttribute(XMLObject.ID) + " but will be ignored, as a data model is provided. Please check the data model for this attribute.");
            }
        }

        boolean integratedProcess = dataModel instanceof IntegratedProcessTraitDataModel; //TODO: can add to interface if that would be better
        if (branchProcessSpecification != null) {
            if (implementationMode == ImplementationMode.LEGACY) {
                throw new XMLParseException("<" + CANONICAL_OU_BRANCH_PROCESSES
                        + "> requires implementation=\"canonical-message\" or implementation=\"canonical\".");
            }
            if (integratedProcess) {
                throw new XMLParseException("<" + CANONICAL_OU_BRANCH_PROCESSES
                        + "> is not yet supported for integrated-process OU trait data.");
            }
        }

        if (delegateProvider == DelegateProvider.OU) {
            if (!integratedProcess) {
                diffusionProcessDelegate = branchProcessSpecification == null
                        ? new OUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels, elasticModel)
                        : new OUDiffusionModelDelegate(
                                treeModel,
                                diffusionModel,
                                optimalTraitsModels,
                                branchProcessSpecification);
            } else {
                diffusionProcessDelegate = new IntegratedOUDiffusionModelDelegate(treeModel, diffusionModel, optimalTraitsModels, elasticModel);
            }
        } else if (delegateProvider == DelegateProvider.DRIFT) {
            diffusionProcessDelegate = new DriftDiffusionModelDelegate(treeModel, diffusionModel, driftModels);
        } else if (delegateProvider == DelegateProvider.HOMOGENOUS) {
            diffusionProcessDelegate = new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        } else {
            throw new RuntimeException("Unrecognized delegate provider.");
        }
        validateCanonicalOUOptions(
                canonicalOUOptions,
                implementationMode,
                implementation,
                delegateProvider);

        // End Parse Evolution Model


        traitName = xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME);
        dataModel.setTipTraitName(getTipTraitName(traitName)); // TODO: not an ideal solution as the trait name could be set differently later

        ConjugateRootTraitPrior rootPrior = ConjugateRootTraitPrior.parseConjugateRootTraitPrior(xo, dataModel.getTraitDimension());

        final boolean allowSingular;

        if (xo.hasAttribute(ALLOW_SINGULAR)) {
            //TODO: check compatibility (there are cases where allowSingular=false guarantees the incorrect likelihood)
            allowSingular = xo.getBooleanAttribute(ALLOW_SINGULAR);
        } else {
            allowSingular = dataModel.getDefaultAllowSingular();
        }
        // End Parse Data Model

        dataModel.addTreeAndRateModel(treeModel, rateTransformation);


        boolean forceCompletelyObserved = xo.getAttribute(FORCE_COMPLETELY_OBSERVED, false);
        final boolean reconstructTraits = xo.getAttribute(RECONSTRUCT_TRAITS, true);
        final boolean initializeLegacyCdi = implementationMode != ImplementationMode.MOMENT_BM;
        // The old ContinuousDataLikelihoodDelegate only ever runs the direct legacy/SafeIntegrator
        // CDI path (see its calculateLikelihood()); it cannot support any of the newer canonical/
        // moment-form backends (its enableCanonicalLikelihood()/enableCdiBMLikelihood()/
        // enableMomentBMLikelihood()/enableCanonicalOULikelihood() all throw
        // UnsupportedOperationException by design). Those backends live only on
        // ContinuousDataLikelihoodDelegateNew, so any non-legacy implementationMode must construct
        // that class instead, or applyImplementationMode's enable*Likelihood() call below always
        // throws.
        //
        // ContinuousDataLikelihoodDelegateNew's own initializeLegacyCdi=true branch is itself
        // deprecated (it unconditionally throws UnsupportedOperationException) since none of its
        // enable*Likelihood() methods need an actual legacy cdi/SafeIntegrator instance - they all
        // build their engine via applyEngineSelection on top of the initializeLegacyCdi=false
        // baseline. So DelegateNew must always be constructed with initializeLegacyCdi=false,
        // regardless of implementationMode; only the old delegate uses the computed flag above.
        final CtsCompatibility delegate = implementationMode == ImplementationMode.LEGACY
                ? new ContinuousDataLikelihoodDelegate(treeModel,
                        diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel,
                        forceCompletelyObserved, allowSingular, initializeLegacyCdi)
                : new ContinuousDataLikelihoodDelegateNew(treeModel,
                        diffusionProcessDelegate, dataModel, rootPrior, rateTransformation, rateModel,
                        forceCompletelyObserved, allowSingular, false);

        if (implementationMode!=ImplementationMode.LEGACY) {
            applyImplementationMode(
                    implementationMode,
                    implementation,
                    delegateProvider,
                    delegate,
                    reconstructTraits,
                    canonicalOUOptions);
        }



        if (dataModel instanceof IntegratedFactorAnalysisLikelihood) {
            if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
                throw new XMLParseException("traitDataLikelihood implementation=\"" + implementation
                        + "\" does not yet support a factor-analysis data model; use implementation=\"legacy\".");
            }
            ((IntegratedFactorAnalysisLikelihood) dataModel).setLikelihoodDelegate((ContinuousDataLikelihoodDelegate) delegate); //TODO
        }

        TreeDataLikelihood treeDataLikelihood = new TreeDataLikelihood(delegate, treeModel,
                rateModel);
        // End Assemble Model

        // Begin Trait Reconstruction Parsing
        if (reconstructTraits) {

//            if (missingIndices != null && missingIndices.size() == 0) {

            ProcessSimulationDelegate simulationDelegate =
                    delegate.getPrecisionType() == PrecisionType.SCALAR ?
                            new ConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                    diffusionModel, dataModel, rootPrior, rateTransformation, delegate) :
                            new MultivariateConditionalOnTipsRealizedDelegate(traitName, treeModel,
                                    diffusionModel, dataModel, rootPrior, rateTransformation, delegate);

            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, simulationDelegate);

            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

            if (dataModel.usesMissingIndices()) {

                ProcessSimulationDelegate fullConditionalDelegate = new TipRealizedValuesViaFullConditionalDelegate(
                        traitName, treeModel, diffusionModel, dataModel, rootPrior, rateTransformation, delegate);

                treeDataLikelihood.addTraits(new ProcessSimulation(treeDataLikelihood, fullConditionalDelegate).getTreeTraits());
            }

            //TODO: remove below (should let ConditionalTraitSimulationHelper figure everything out)
            int[] partitionDimensions = dataModel.getPartitionDimensions();
            if (partitionDimensions.length > 1) {
                PartitionedTreeTraitProvider partitionedProvider =
                        new PartitionedTreeTraitProvider(treeDataLikelihood.getTreeTraits(), partitionDimensions);
                treeDataLikelihood.addTraits(partitionedProvider.getTreeTraits());
            }

        }
        // End Trait Reconstruction Parsing

        delegate.setExtensionHelper();

        return treeDataLikelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of trait data on a tree given a diffusion model.";
    }

    public Class getReturnType() {
        return TreeDataLikelihood.class;
    }

    private static void applyImplementationMode(final ImplementationMode implementationMode,
                                                final String implementation,
                                                final DelegateProvider delegateProvider,
                                                final CtsCompatibility delegate,
                                                final boolean reconstructTraits,
                                                final CanonicalOUOptions canonicalOUOptions) throws XMLParseException {
        if (implementationMode.requiresLikelihoodOnly() && reconstructTraits) {
            throw new XMLParseException(
                    "traitDataLikelihood implementation=\"" + implementation
                            + "\" does not yet support reconstructTraits=\"true\".");
        }

        try {
            switch (implementationMode) {
                case LEGACY:
                    return;
                case CANONICAL_COMPATIBILITY:
                    delegate.enableCanonicalLikelihood(canonicalOUOptions);
                    return;
                case CDI_MOMENT_BM:
                    if (delegateProvider != DelegateProvider.HOMOGENOUS) {
                        throw new XMLParseException(
                                "traitDataLikelihood implementation=\"" + implementation
                                        + "\" is the CDI/moment-form BM production backend and supports only homogeneous BM. "
                                        + "Use implementation=\"legacy\" for legacy CDI OU/drift models or "
                                        + "implementation=\"canonical-message\" for canonical-message OU.");
                    }
                    delegate.enableCdiBMLikelihood();
                    return;
                case MOMENT_BM:
                    if (delegateProvider != DelegateProvider.HOMOGENOUS) {
                        throw new XMLParseException(
                                "traitDataLikelihood implementation=\"" + implementation
                                        + "\" is the native moment-form BM backend and supports only homogeneous BM. "
                                        + "Use implementation=\"legacy\" for legacy CDI OU/drift models or "
                                        + "implementation=\"canonical-message\" for canonical-message OU.");
                    }
                    delegate.enableMomentBMLikelihood();
                    return;
                case CANONICAL_MESSAGE_OU:
                    if (delegateProvider != DelegateProvider.OU) {
                        throw new XMLParseException(
                                "traitDataLikelihood implementation=\"" + implementation
                                        + "\" is the canonical-message OU backend and supports only OU models. "
                                        + "Use implementation=\"cdi\" or \"cdi-moment\" for production BM.");
                    }
                    delegate.enableCanonicalOULikelihood(canonicalOUOptions);
                    return;
                case NATURAL_CANONICAL_BM:
                    throw new XMLParseException(
                            "traitDataLikelihood implementation=\"" + implementation
                                    + "\" names the experimental natural-parameter BM backend. "
                                    + "It is retained for reference only and is not a production XML backend; "
                                    + "use implementation=\"cdi\" or \"cdi-moment\" for BM.");
                default:
                    throw new XMLParseException(
                            "Unknown traitDataLikelihood implementation mode for \"" + implementation + "\".");
            }
        } catch (UnsupportedOperationException ex) {
            throw new XMLParseException(
                    "traitDataLikelihood implementation=\"" + implementation
                            + "\" is not supported here: " + ex.getMessage());
        }
    }

    private static String quotedImplementationList() {
        return "\"" + IMPLEMENTATION_LEGACY + "\", "
                + "\"" + IMPLEMENTATION_CANONICAL + "\", "
                + "\"" + IMPLEMENTATION_CDI + "\", "
                + "\"" + IMPLEMENTATION_CDI_MOMENT + "\", "
                + "\"" + IMPLEMENTATION_BM_CDI + "\", "
                + "\"" + IMPLEMENTATION_BM_CDI_MOMENT + "\", "
                + "\"" + IMPLEMENTATION_BM_MOMENT + "\", "
                + "\"" + IMPLEMENTATION_NATIVE_MOMENT_BM + "\", "
                + "\"" + IMPLEMENTATION_CANONICAL_MESSAGE + "\", "
                + "\"" + IMPLEMENTATION_OU_CANONICAL + "\", "
                + "\"" + IMPLEMENTATION_OU_CANONICAL_MESSAGE + "\", "
                + "\"" + IMPLEMENTATION_NATURAL_CANONICAL_BM + "\", or "
                + "\"" + IMPLEMENTATION_BM_NATURAL_CANONICAL + "\"";
    }

    private static CanonicalOUOptions parseCanonicalOUOptions(final XMLObject xo)
            throws XMLParseException {
        final String coordinates = xo.getAttribute(
                CANONICAL_OU_COORDINATES,
                CanonicalOUOptions.COORDINATES_AUTO);
        try {
            return CanonicalOUOptions.parseCoordinates(coordinates);
        } catch (IllegalArgumentException ex) {
            throw new XMLParseException(ex.getMessage());
        }
    }

    private static void validateCanonicalOUOptions(final CanonicalOUOptions canonicalOUOptions,
                                                   final ImplementationMode implementationMode,
                                                   final String implementation,
                                                   final DelegateProvider delegateProvider)
            throws XMLParseException {
        if (!canonicalOUOptions.isExplicit()) {
            return;
        }
        final boolean selectsCanonicalOU =
                delegateProvider == DelegateProvider.OU
                        && (implementationMode == ImplementationMode.CANONICAL_COMPATIBILITY
                        || implementationMode == ImplementationMode.CANONICAL_MESSAGE_OU);
        if (!selectsCanonicalOU) {
            throw new XMLParseException(
                    CANONICAL_OU_COORDINATES + "=\"" + canonicalOUOptions.toXmlString()
                            + "\" is only valid for canonical OU traitDataLikelihood backends; "
                            + "implementation=\"" + implementation + "\" selected "
                            + delegateProvider.name().toLowerCase(Locale.ROOT) + ".");
        }
    }

    private static CanonicalOUBranchProcessSpecification parseCanonicalOUBranchProcesses(
            final XMLObject xo) throws XMLParseException {
        if (!xo.hasChildNamed(CANONICAL_OU_BRANCH_PROCESSES)) {
            return null;
        }

        final XMLObject cxo = xo.getChild(CANONICAL_OU_BRANCH_PROCESSES);
        if (!cxo.hasChildNamed(BRANCH_PROCESS_INDEX)) {
            throw new XMLParseException("<" + CANONICAL_OU_BRANCH_PROCESSES
                    + "> requires a <" + BRANCH_PROCESS_INDEX + "> child.");
        }

        final XMLObject branchProcessIndexXO = cxo.getChild(BRANCH_PROCESS_INDEX);
        final Parameter branchProcessIndex =
                (Parameter) branchProcessIndexXO.getChild(Parameter.class);
        if (branchProcessIndex == null) {
            throw new XMLParseException("<" + BRANCH_PROCESS_INDEX
                    + "> must contain a branch-process index parameter.");
        }

        final List<MultivariateElasticModel> elasticModels =
                new ArrayList<MultivariateElasticModel>();
        for (int i = 0; i < cxo.getChildCount(); ++i) {
            final Object child = cxo.getChild(i);
            if (!(child instanceof XMLObject)) {
                continue;
            }
            final XMLObject processXO = (XMLObject) child;
            if (!PROCESS.equals(processXO.getName())) {
                continue;
            }
            elasticModels.add(parseBranchProcessElasticModel(processXO));
        }
        if (elasticModels.isEmpty()) {
            throw new XMLParseException("<" + CANONICAL_OU_BRANCH_PROCESSES
                    + "> requires at least one <" + PROCESS + "> child.");
        }
        return new CanonicalOUBranchProcessSpecification(elasticModels, branchProcessIndex);
    }

    private static MultivariateElasticModel parseBranchProcessElasticModel(
            final XMLObject processXO) throws XMLParseException {
        if (!processXO.hasChildNamed(STRENGTH_OF_SELECTION_MATRIX)) {
            throw new XMLParseException("<" + PROCESS + "> inside <"
                    + CANONICAL_OU_BRANCH_PROCESSES + "> requires <"
                    + STRENGTH_OF_SELECTION_MATRIX + ">.");
        }
        final XMLObject selectionXO = processXO.getChild(STRENGTH_OF_SELECTION_MATRIX);
        final MatrixParameterInterface strengthOfSelectionMatrix =
                (MatrixParameterInterface) selectionXO.getChild(MatrixParameterInterface.class);
        if (strengthOfSelectionMatrix == null) {
            throw new XMLParseException("<" + STRENGTH_OF_SELECTION_MATRIX
                    + "> inside <" + PROCESS + "> must contain a matrix parameter.");
        }
        return new MultivariateElasticModel(strengthOfSelectionMatrix);
    }

    private static void validateBranchProcessSelectionCharts(
            final XMLObject xo,
            final CanonicalOUBranchProcessSpecification specification) throws XMLParseException {
        for (int process = 0; process < specification.getProcessCount(); ++process) {
            OUSelectionChartParserHelper.validateSelectionChart(
                    xo,
                    specification.getElasticModel(process).getStrengthOfSelectionMatrixParameter(),
                    CONTINUOUS_DATA_LIKELIHOOD + " " + CANONICAL_OU_BRANCH_PROCESSES
                            + " process " + process);
        }
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(MultivariateDiffusionModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(CONJUGATE_ROOT_PRIOR, ConjugateRootTraitPrior.rules),
            new XORRule(new XMLSyntaxRule[]{
                    new AndRule(ContinuousTraitDataModelParser.rules),
                    new ElementRule(ContinuousTraitPartialsProvider.class)
            }),
            new ElementRule(DRIFT_MODELS, new XMLSyntaxRule[]{
                    new ElementRule(BranchRateModel.class, 1, Integer.MAX_VALUE),
            }, true),
            new ElementRule(OPTIMAL_TRAITS, new XMLSyntaxRule[]{
                    new ElementRule(BranchRateModel.class, 1, Integer.MAX_VALUE),
            }, true),
            new ElementRule(STRENGTH_OF_SELECTION_MATRIX,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}, true),
            new ElementRule(CANONICAL_OU_BRANCH_PROCESSES, new XMLSyntaxRule[]{
                    new ElementRule(BRANCH_PROCESS_INDEX,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(PROCESS, new XMLSyntaxRule[]{
                            new ElementRule(STRENGTH_OF_SELECTION_MATRIX,
                                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)})
                    }, 1, Integer.MAX_VALUE)
            }, true),
            AttributeRule.newBooleanRule(SCALE_BY_TIME, true),
            AttributeRule.newBooleanRule(USE_TREE_LENGTH, true),
            AttributeRule.newBooleanRule(RECIPROCAL_RATES, true),
            AttributeRule.newBooleanRule(RECONSTRUCT_TRAITS, true),
            AttributeRule.newBooleanRule(ALLOW_SINGULAR, true),
            AttributeRule.newBooleanRule(FORCE_FULL_PRECISION, true),
            AttributeRule.newBooleanRule(FORCE_COMPLETELY_MISSING, true),
            AttributeRule.newBooleanRule(FORCE_DRIFT, true),
            AttributeRule.newBooleanRule(FORCE_OU, true),
            AttributeRule.newBooleanRule(FORCE_COMPLETELY_OBSERVED, true),
            AttributeRule.newBooleanRule(REPORT_MULTIVARIATE, true),
            AttributeRule.newBooleanRule(INTEGRATE_INTERNAL_TRAITS, true),
            new StringAttributeRule(OUSelectionChartParserHelper.SELECTION_CHART,
                    "Selection-matrix chart for OU models. Orthogonal block is the default; dense must be explicit.",
                    OUSelectionChartParserHelper.ALLOWED_SELECTION_CHARTS, true),
            AttributeRule.newStringRule(IMPLEMENTATION, true),
            new StringAttributeRule(CANONICAL_OU_COORDINATES,
                    "Coordinate system for canonical OU messages. Use homogeneousBlockBasis "
                            + "for homogeneous block-diagonal OU in the selection basis.",
                    CanonicalOUOptions.ALLOWED_COORDINATES,
                    true),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            TreeTraitParserUtilities.jitterRules(true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    //********************************************************************
    // BeautiModelIDProvider interface
    //********************************************************************

    private static final String CONTINUOUS_DATA_LIKELIHOOD_DEFAULT_ID = "traitLikelihood";

    public String getParserTag() {
        return CONTINUOUS_DATA_LIKELIHOOD;
    }

    public String getId(String modelName) {
        return modelName + "." + CONTINUOUS_DATA_LIKELIHOOD_DEFAULT_ID;
    }
}
