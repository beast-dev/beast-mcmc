/*
 * BeagleDataLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.basta;

/**
 * BeagleDataLikelihoodDelegate
 *
 * A DataLikelihoodDelegate that uses BEAGLE
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

import beagle.*;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.EvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.HomogenousSubstitutionModelDelegate;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static dr.evomodel.treedatalikelihood.BeagleFunctionality.*;

@Deprecated
public class BeagleMatrixExponentiationDelegate extends AbstractModel implements Citable {

    // This property is a comma-delimited list of resource numbers (0 == CPU) to
    // allocate each BEAGLE instance to. If less than the number of instances then
    // will wrap around.
    private static final String RESOURCE_AUTO_PROPERTY = "beagle.resource.auto";
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    private static final String SCALING_PROPERTY = "beagle.scaling";
    private static final String RESCALE_FREQUENCY_PROPERTY = "beagle.rescale";
    private static final String DELAY_SCALING_PROPERTY = "beagle.delay.scaling";
    private static final String EXTRA_BUFFER_COUNT_PROPERTY = "beagle.extra.buffer.count";
    private static final String FORCE_VECTORIZATION = "beagle.force.vectorization";
    private static final String THREAD_COUNT = "beagle.thread.count";

    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private static List<String> scalingOrder = null;
    private static List<Integer> extraBufferOrder = null;

    private static final boolean DEBUG = false;

    /**
     *
     */
    public BeagleMatrixExponentiationDelegate(SubstitutionModel substitutionModel, int branchCount) {

        super("BeagleMatrixExponentiationDelegate");
        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing BEAGLE Matrix Exponentiation Delegate");

        this.substitutionModel = substitutionModel;

        this.dataType = this.substitutionModel.getDataType();
        stateCount = dataType.getStateCount();

        addModel(this.substitutionModel);

        this.branchCount = branchCount;

        branchUpdates = new boolean[this.branchCount];
        storedBranchUpdates = new boolean[this.branchCount];
        branchLengths = new double[this.branchCount];
        storedBranchLengths = new double[this.branchCount];

        evolutionaryProcessDelegate = new HomogenousSubstitutionModelDelegate(
                null, new HomogeneousBranchModel(substitutionModel));

        // Attempt to get the resource order from the System Property
        if (resourceOrder == null) {
            resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (preferredOrder == null) {
            preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (requiredOrder == null) {
            requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }
        if (extraBufferOrder == null) {
            extraBufferOrder = parseSystemPropertyIntegerArray(EXTRA_BUFFER_COUNT_PROPERTY);
        }

        int[] resourceList = null;
        long preferenceFlags = 0;
        long requirementFlags = 0;

        if (resourceOrder.size() > 0) {
            // added the zero on the end so that a CPU is selected if requested resource fails
            resourceList = new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
            if (resourceList[0] > 0) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }
        }

        if (preferredOrder.size() > 0) {
            preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
        }

        if (requiredOrder.size() > 0) {
            requirementFlags = requiredOrder.get(instanceCount % requiredOrder.size());
        }

        boolean forceVectorization = false;
        String vectorizationString = System.getProperty(FORCE_VECTORIZATION);
        if (vectorizationString != null) {
            forceVectorization = true;
        }

        int threadCount = -1;
        String tc = System.getProperty(THREAD_COUNT);
        if (tc != null) {
            threadCount = Integer.parseInt(tc);
        }

        if (threadCount == 0 || threadCount == 1) {
            preferenceFlags &= ~BeagleFlag.THREADING_CPP.getMask();
            preferenceFlags |= BeagleFlag.THREADING_NONE.getMask();
        } else {
            preferenceFlags &= ~BeagleFlag.THREADING_NONE.getMask();
            preferenceFlags |= BeagleFlag.THREADING_CPP.getMask();
        }

        if (BeagleFlag.VECTOR_SSE.isSet(preferenceFlags) && (stateCount != 4)
                && !forceVectorization && !IS_ODD_STATE_SSE_FIXED()
                ) {
            // @todo SSE doesn't seem to work for larger state spaces so for now we override the
            // SSE option.
            preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
            preferenceFlags |= BeagleFlag.VECTOR_NONE.getMask();
        }

        if (!BeagleFlag.PRECISION_SINGLE.isSet(preferenceFlags)) {
            // if single precision not explicitly set then prefer double
            preferenceFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
        }

        if (evolutionaryProcessDelegate.canReturnComplexDiagonalization()) {
            requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();
        }

        if ((resourceList == null &&
                (BeagleFlag.PROCESSOR_GPU.isSet(preferenceFlags) ||
                        BeagleFlag.FRAMEWORK_CUDA.isSet(preferenceFlags) ||
                        BeagleFlag.FRAMEWORK_OPENCL.isSet(preferenceFlags)))
                ||
                (resourceList != null && resourceList[0] > 0)) {
            // non-CPU implementations don't have SSE so remove default preference for SSE
            // when using non-CPU preferences or prioritising non-CPU resource
            preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
            preferenceFlags &= ~BeagleFlag.THREADING_CPP.getMask();
        }

        beagle = BeagleFactory.loadBeagleInstance(
                0,
                0,
                0,
                stateCount,
                0,
                evolutionaryProcessDelegate.getEigenBufferCount(),
                evolutionaryProcessDelegate.getMatrixBufferCount(),
                0,
                0, // Always allocate; they may become necessary
                resourceList,
                preferenceFlags,
                requirementFlags
        );

        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;

        if (instanceDetails != null) {
            resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
            if (resourceDetails != null) {
                StringBuilder sb = new StringBuilder("  Using BEAGLE resource ");
                sb.append(resourceDetails.getNumber()).append(": ");
                sb.append(resourceDetails.getName()).append("\n");
                if (resourceDetails.getDescription() != null) {
                    String[] description = resourceDetails.getDescription().split("\\|");
                    for (String desc : description) {
                        if (desc.trim().length() > 0) {
                            sb.append("    ").append(desc.trim()).append("\n");
                        }
                    }
                }
                sb.append("    with instance flags: ").append(instanceDetails.toString());
                logger.info(sb.toString());
            } else {
                logger.info("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
            }
        } else {
            logger.info("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
        }

        if (IS_THREAD_COUNT_COMPATIBLE() && threadCount > 1) {
            beagle.setCPUThreadCount(threadCount);
        }

        updateSubstitutionModel = true;


        instanceCount++;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public void calculateTransitionProbabilityMatrices() {

        //recomputeScaleFactors = false;
        if (DEBUG) {
//            System.out.println("Partition: " + this.getModelName());
        }

        int branchUpdateCount = 0;
        int[] branchUpdateIndices = new int[branchCount];
        double[] branchLengths = new double[branchCount];
        for (int i = 0; i < branchCount; i++) {
            if (this.branchUpdates[i]) {
                branchUpdateIndices[branchUpdateCount] = i;
                branchLengths[branchUpdateCount] = branchLengths[i];
                branchUpdateCount++;
            }
        }

        if (updateSubstitutionModel) { // TODO More efficient to update only the substitution model that changed, instead of all
            evolutionaryProcessDelegate.updateSubstitutionModels(beagle, true);
        }

        if (branchUpdateCount > 0) {
            evolutionaryProcessDelegate.updateTransitionMatrices(
                    beagle,
                    branchUpdateIndices,
                    branchLengths,
                    branchUpdateCount,
                    true);
        }

        if (DEBUG) {
//            System.out.println("useScaleFactors=" + useScaleFactors + " recomputeScaleFactors=" + recomputeScaleFactors + " (" + getId() + ")");
        }

        updateSubstitutionModel = false;
    }

    public void makeDirty() {
        updateSubstitutionModel = true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == substitutionModel) {
            updateSubstitutionModel = true;
        }



        // Tell listeners to update
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    /**
     * Stores the additional state other than model components
     */
    @Override
    public void storeState() {
        evolutionaryProcessDelegate.storeState();

        System.arraycopy(branchUpdates, 0, storedBranchUpdates, 0, branchCount);
        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchCount);
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restoreState() {
        evolutionaryProcessDelegate.restoreState();

        boolean[] tmp1 = branchUpdates;
        branchUpdates = storedBranchUpdates;
        storedBranchUpdates = tmp1;

        double[] tmp2 = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp2;

    }

    @Override
    protected void acceptState() {
    }

    // **************************************************************
    // INSTANCE CITABLE
    // **************************************************************

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using BEAGLE likelihood calculation library";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.AYRES_2019_BEAGLE);
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private final int branchCount;

    private boolean[] branchUpdates;
    private double[] branchLengths;

    private boolean[] storedBranchUpdates;
    private double[] storedBranchLengths;


    /**
     * the data type
     */
    private final DataType dataType;

    /**
     * the number of states in the data
     */
    private final int stateCount;

    /**
     * the substitution model
     */
    private final SubstitutionModel substitutionModel;

    /**
     * A delegate to handle substitution models on branches
     */
    private final EvolutionaryProcessDelegate evolutionaryProcessDelegate;

    /**
     * the BEAGLE library instance
     */
    private final Beagle beagle;

    /**
     * Flag to specify that the substitution model has changed
     */
    private boolean updateSubstitutionModel;

}
