/*
 * PartitionedGradients.java
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

package dr.evomodel.treedatalikelihood;

import dr.inference.hmc.CompoundParameterGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * Gradients of a likelihood split into partitions, one TreeDataLikelihood each (a multi-partition
 * {@code <treeDataLikelihood>} parses to a CompoundLikelihood of them on the CPU):
 * <ul>
 * <li>a parameter shared by several partitions gets the sum of its per-partition gradients ({@link JointGradient});</li>
 * <li>a CompoundParameter that no partition uses whole, such as allNus, is differentiated block by block, each block
 * summed over the partitions that use it ({@link CompoundParameterGradient}).</li>
 * </ul>
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class PartitionedGradients {

    private PartitionedGradients() {
    }

    /**
     * Builds the gradient of one partition wrt a parameter it uses.
     */
    public interface Factory {

        boolean appliesTo(TreeDataLikelihood likelihood, Parameter parameter);

        GradientWrtParameterProvider create(TreeDataLikelihood likelihood, Parameter parameter);
    }

    /**
     * @return the likelihood itself, or the members of a CompoundLikelihood, which must all be TreeDataLikelihoods
     */
    public static List<TreeDataLikelihood> treeDataLikelihoods(Likelihood likelihood) {
        List<TreeDataLikelihood> list = new ArrayList<>();
        if (likelihood instanceof TreeDataLikelihood) {
            list.add((TreeDataLikelihood) likelihood);
        } else if (likelihood instanceof CompoundLikelihood) {
            for (Likelihood member : ((CompoundLikelihood) likelihood).getLikelihoods()) {
                if (!(member instanceof TreeDataLikelihood)) {
                    throw new IllegalArgumentException(likelihood.getId() + " contains " + member.getId()
                            + ", which is not a treeDataLikelihood");
                }
                list.add((TreeDataLikelihood) member);
            }
        } else {
            throw new IllegalArgumentException(likelihood.getId() + " is neither a treeDataLikelihood nor a"
                    + " compound of them");
        }
        return list;
    }

    /**
     * @param simplexSum the fixed sum of a compound parameter differentiated block by block, or NaN
     */
    public static GradientWrtParameterProvider create(List<TreeDataLikelihood> likelihoods, Parameter parameter,
                                                      Factory factory, double simplexSum, Double tolerance) {

        if (appliesToAny(likelihoods, parameter, factory)) {
            GradientWrtParameterProvider whole = sum(likelihoods, parameter, factory);
            if (!Double.isNaN(simplexSum) && Double.compare(simplexSum, whole.getSimplexSum()) != 0) {
                throw new IllegalArgumentException("simplexSum is only for a compound parameter whose blocks belong to"
                        + " different partitions, not " + parameter.getId());
            }
            return whole;
        }

        if (parameter instanceof CompoundParameter) {
            CompoundParameter compound = (CompoundParameter) parameter;
            List<GradientWrtParameterProvider> blocks = new ArrayList<>();
            for (int i = 0; i < compound.getParameterCount(); ++i) {
                blocks.add(sum(likelihoods, compound.getParameter(i), factory));
            }
            return new CompoundParameterGradient(compound, blocks, simplexSum, tolerance);
        }

        throw new IllegalArgumentException(parameter.getId() + " is not used by any partition");
    }

    private static boolean appliesToAny(List<TreeDataLikelihood> likelihoods, Parameter parameter, Factory factory) {
        for (TreeDataLikelihood likelihood : likelihoods) {
            if (factory.appliesTo(likelihood, parameter)) {
                return true;
            }
        }
        return false;
    }

    private static GradientWrtParameterProvider sum(List<TreeDataLikelihood> likelihoods, Parameter parameter,
                                                    Factory factory) {
        List<GradientWrtParameterProvider> providers = new ArrayList<>();
        for (TreeDataLikelihood likelihood : likelihoods) {
            if (factory.appliesTo(likelihood, parameter)) {
                providers.add(factory.create(likelihood, parameter));
            }
        }
        if (providers.isEmpty()) {
            throw new IllegalArgumentException(parameter.getId() + " is not used by any partition");
        }
        return providers.size() == 1 ? providers.get(0) : new JointGradient(providers);
    }
}
