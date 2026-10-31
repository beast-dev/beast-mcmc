/*
 * ThreadingScheme.java
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

package dr.app.beast;

import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.mcmc.MCMC;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Decides how each chain uses its share of the processors (its thread budget):
 *
 * <ul>
 * <li>serial: the chain's likelihoods are computed one after another and every BEAGLE CPU
 * instance gets the whole budget.</li>
 * <li>split: the likelihoods are computed in parallel and the budget is divided between the
 * BEAGLE CPU instances.</li>
 * <li>auto: split if the chain has at least two BEAGLE GPU instances, or at least two BEAGLE
 * CPU instances that each have enough work to be worth a thread of their own (and a budget of
 * at least two threads to share); serial otherwise.</li>
 * </ul>
 *
 * Benchmarks (see documentation/multi-chain.md) found computing likelihoods in parallel worth
 * it only for several partitions that each have substantial data, and costly for small data;
 * giving all of the budget to BEAGLE and all of it to the parallel likelihoods at the same time
 * (the old behaviour) was never the fastest.
 *
 * The scheme is applied to each chain as it starts, when its model (and so its BEAGLE
 * instances) exist. A thread count given explicitly (-threads or -beagle_threads) is kept.
 *
 * @author Andrew Rambaut
 */
public class ThreadingScheme implements MCMC.ChainPreparer {

    /**
     * Below this amount of work (tips x patterns x (states / 4)^2) a BEAGLE CPU instance is not
     * worth a thread of its own: a 165 taxon, 90 pattern partition (15,000) was much slower
     * computed in parallel; 1441 x 200 (288,000) and 899 x 430 (386,000) were faster.
     */
    public static final double MIN_WORK_PER_INSTANCE = 100000;

    public enum Type {
        SERIAL("serial"),
        SPLIT("split"),
        AUTO("auto");

        Type(String name) {
            this.name = name;
        }

        public static Type parse(String name) {
            for (Type type : values()) {
                if (type.name.equalsIgnoreCase(name)) {
                    return type;
                }
            }
            throw new IllegalArgumentException("Unknown threading scheme '" + name + "' (use serial, split or auto)");
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }

    /**
     * @param type                  the scheme
     * @param threadBudget          the number of threads each chain may use
     * @param fixedLikelihoodThreads true if the likelihood thread pools were set with -threads
     * @param fixedBeagleThreads    true if the BEAGLE CPU thread count was set explicitly
     */
    public ThreadingScheme(Type type, int threadBudget, boolean fixedLikelihoodThreads, boolean fixedBeagleThreads) {
        this.type = type;
        this.threadBudget = Math.max(1, threadBudget);
        this.fixedLikelihoodThreads = fixedLikelihoodThreads;
        this.fixedBeagleThreads = fixedBeagleThreads;
    }

    @Override
    public void prepare(MCMC mcmc) {
        Likelihood likelihood = mcmc.getLikelihood();

        Set<Likelihood> all = new LinkedHashSet<Likelihood>();
        all.add(likelihood);
        all.addAll(likelihood.getLikelihoodSet());

        List<BeagleDataLikelihoodDelegate> cpuInstances = new ArrayList<BeagleDataLikelihoodDelegate>();
        int gpuInstanceCount = 0;
        int heavyCPUInstanceCount = 0;
        List<CompoundLikelihood> compounds = new ArrayList<CompoundLikelihood>();

        for (Likelihood l : all) {
            if (l instanceof CompoundLikelihood) {
                compounds.add((CompoundLikelihood) l);
            } else if (l instanceof TreeDataLikelihood) {
                DataLikelihoodDelegate delegate = ((TreeDataLikelihood) l).getDataLikelihoodDelegate();
                if (delegate instanceof BeagleDataLikelihoodDelegate) {
                    BeagleDataLikelihoodDelegate beagle = (BeagleDataLikelihoodDelegate) delegate;
                    if (beagle.isCPUInstance()) {
                        cpuInstances.add(beagle);
                        if (getWork(beagle) >= MIN_WORK_PER_INSTANCE) {
                            heavyCPUInstanceCount++;
                        }
                    } else {
                        gpuInstanceCount++;
                    }
                } else if (delegate instanceof MultiPartitionDataLikelihoodDelegate) {
                    // the multi-partition delegate only runs on GPUs
                    gpuInstanceCount++;
                }
            }
        }

        Type chosen = type;
        if (type == Type.AUTO) {
            boolean split = gpuInstanceCount >= 2 || (heavyCPUInstanceCount >= 2 && threadBudget >= 2);
            chosen = split ? Type.SPLIT : Type.SERIAL;
        }

        int beagleThreads = (chosen == Type.SPLIT ?
                Math.max(1, threadBudget / Math.max(1, cpuInstances.size())) : threadBudget);

        int likelihoodThreads = 0;
        if (!fixedLikelihoodThreads) {
            // With split, a compound likelihood that directly contains several BEAGLE
            // likelihoods (or compounds of them) computes them in parallel, with a thread each.
            // The data likelihoods are usually in a compound of their own inside the posterior.
            for (CompoundLikelihood compound : compounds) {
                int threads = 0;
                if (chosen == Type.SPLIT) {
                    for (Likelihood child : compound.getLikelihoods()) {
                        if (containsBeagle(child)) {
                            threads++;
                        }
                    }
                }
                compound.setThreadCount(threads > 1 ? threads : 0);
                likelihoodThreads = Math.max(likelihoodThreads, threads > 1 ? threads : 0);
            }
        }
        int changed = 0;
        if (!fixedBeagleThreads) {
            for (BeagleDataLikelihoodDelegate beagle : cpuInstances) {
                if (beagle.setCPUThreadCount(beagleThreads)) {
                    changed++;
                }
            }
        }

        if (REPORTED.compareAndSet(false, true)) {
            StringBuilder sb = new StringBuilder("Threading scheme: ").append(chosen);
            if (type == Type.AUTO) {
                sb.append(" (auto)");
            }
            sb.append(", ").append(threadBudget).append(" thread").append(threadBudget > 1 ? "s" : "").append(" per chain\n");
            if (fixedLikelihoodThreads) {
                sb.append("  likelihood threads as set by -threads\n");
            } else if (likelihoodThreads > 1) {
                sb.append("  likelihoods computed in parallel (").append(likelihoodThreads).append(" threads)\n");
            } else {
                sb.append("  likelihoods computed one after another\n");
            }
            sb.append("  BEAGLE: ").append(cpuInstances.size()).append(" CPU instance").append(cpuInstances.size() != 1 ? "s" : "");
            if (fixedBeagleThreads) {
                sb.append(" (threads as set by -beagle_threads)");
            } else if (!cpuInstances.isEmpty()) {
                sb.append(" with ").append(beagleThreads).append(" thread").append(beagleThreads > 1 ? "s" : "")
                        .append(cpuInstances.size() > 1 ? " each" : "");
                if (changed < cpuInstances.size()) {
                    sb.append(" (").append(cpuInstances.size() - changed).append(" without threading)");
                }
            }
            sb.append(", ").append(gpuInstanceCount).append(" GPU instance").append(gpuInstanceCount != 1 ? "s" : "");
            System.out.println(sb);
        }
    }

    /**
     * @return true if the likelihood is computed by BEAGLE or contains one that is
     */
    private static boolean containsBeagle(Likelihood likelihood) {
        if (isBeagle(likelihood)) {
            return true;
        }
        for (Likelihood l : likelihood.getLikelihoodSet()) {
            if (isBeagle(l)) {
                return true;
            }
        }
        return false;
    }

    private static boolean isBeagle(Likelihood likelihood) {
        if (likelihood instanceof TreeDataLikelihood) {
            DataLikelihoodDelegate delegate = ((TreeDataLikelihood) likelihood).getDataLikelihoodDelegate();
            return delegate instanceof BeagleDataLikelihoodDelegate || delegate instanceof MultiPartitionDataLikelihoodDelegate;
        }
        return false;
    }

    /**
     * @return a measure of the work of one likelihood evaluation, relative to nucleotide data
     */
    static double getWork(BeagleDataLikelihoodDelegate beagle) {
        double states = beagle.getInstanceStateCount() / 4.0;
        return (double) beagle.getTipCount() * beagle.getInstancePatternCount() * states * states;
    }

    public Type getType() {
        return type;
    }

    public int getThreadBudget() {
        return threadBudget;
    }

    private static final AtomicBoolean REPORTED = new AtomicBoolean(false);

    private final Type type;
    private final int threadBudget;
    private final boolean fixedLikelihoodThreads;
    private final boolean fixedBeagleThreads;
}
