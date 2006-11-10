/*
 * Genome.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.simulator;

import dr.simulator.fitnessfunctions.FitnessFunction;

import java.util.*;


/**
 * @author Andrew Rambaut
 */
public class Genome {

	public Genome() { }

	public Genome(byte[] sequence) {

        List mutationList = new ArrayList();
        if (masterSequence == null) {
            throw new IllegalArgumentException("The master sequence has not been set");
        } else {

            if (sequence.length != masterSequence.length) {
                throw new IllegalArgumentException("Initializing sequence length is different to the master sequence");
            }

            for (int i = 0; i < sequence.length; i++) {
                if (sequence[i] != masterSequence[i]) {
                    mutationList.add(Mutation.getMutation(i, sequence[i]));
                }
            }
        }
        this.mutations = new Mutation[mutationList.size()];
        mutationList.toArray(mutations);
	}

    public void duplicate(Genome source) {
        mutations = source.mutations;
        mutationCount = source.mutationCount;
        logFitness = source.logFitness;
    }

    /**
     * @return the length of the genome
     */
	public int getLength() {
		return masterSequence.length;
	}

    /**
     * Gets the state at a given position. If no mutation exists then the masterSequence state is returned.
     * @param position the position in the sequence
     * @return the state
     */
    public byte getState(int position) {
        Mutation mutation = getMutation(position);
        if (mutation != null) {
            return mutation.state;
        }
		return masterSequence[position];
    }

    public int getMutationCount() {
        return mutationCount;
    }

    public Mutation getMutation(int position) {

        if (mutationCount == 0) return null;

        int i0 = 0;
        int i1 = mutationCount;


        do {
            int j = (i1 + i0) / 2;
            if (position < mutations[j].position) {
                i1 = j;
            } else if (position < mutations[j].position) {
                i0 = j;
            } else {
                return mutations[j];
            }
        } while (i0 < i1);

        return null;
    }

    public Mutation[] getMutations() {

        return mutations;
    }

	public void clearMutations() {
	    mutations = null;
	}

    public void applyMutations(Mutation[] newMutations, FitnessFunction fitnessFunction) {

        Mutation[] oldMutations = mutations;
        int oldMutationCount = mutationCount;

        // this might well be more than is required...
        int capacity = oldMutationCount + newMutations.length;

        mutations = new Mutation[capacity];

        int i = 0;
        int j = 0;
        mutationCount = 0;

        int position;
        byte oldState;
        byte newState;

        while (i < oldMutationCount && j < newMutations.length) {
            if (oldMutations[i].position < newMutations[j].position) {
                position = oldMutations[i].position;
                oldState = oldMutations[i].state;
                newState = oldMutations[i].state;
                i++;
            } else if (oldMutations[i].position > newMutations[j].position) {
                position = newMutations[j].position;
                oldState = masterSequence[position];
                newState = newMutations[j].state;
                j++;

                if (newState != oldState) {
                    logFitness -= fitnessFunction.getLogFitness(position, oldState);
                    logFitness += fitnessFunction.getLogFitness(position, newState);
                }
            } else { // they are the same - use the new mutation
                position = newMutations[j].position;
                oldState = oldMutations[i].state;
                newState = newMutations[j].state;
                i++;
                j++;

                if (newState != oldState) {
                    logFitness -= fitnessFunction.getLogFitness(position, oldState);
                    logFitness += fitnessFunction.getLogFitness(position, newState);
                }
            }

            if (newState != masterSequence[position]) {
                // if the new mutation is not a reversion to the masterSequence then add it
                mutations[mutationCount] = Mutation.getMutation(position, newState);
                mutationCount++;
            }
        }

        while (i < oldMutationCount) {
            // add the remaining old mutations (these don't affect the logFitness
            mutations[mutationCount] = oldMutations[i];
            i++;
            mutationCount++;
        }

        while (j < newMutations.length) {
            position = newMutations[j].position;
            newState = newMutations[j].state;
            j++;
            if (newState != masterSequence[position]) {
                // if the new mutation is not a reversion to the masterSequence then add it

                logFitness -= fitnessFunction.getLogFitness(position, masterSequence[position]);
                logFitness += fitnessFunction.getLogFitness(position, newState);

                mutations[mutationCount] = Mutation.getMutation(position, newState);
                mutationCount++;
            }
        }

    }

    public void changeMasterSequence(byte[] newMasterSequence) {

        List mutationList = new ArrayList();
        for (int position = 0; position < masterSequence.length; position++) {
            byte state = getState(position);
            if (state != newMasterSequence[position]) {
                // if the current state differs from the new masterSequence
                // then add a mutation
                mutationList.add(Mutation.getMutation(position, state));
            }
        }
        this.mutations = new Mutation[mutationList.size()];
        mutationList.toArray(mutations);
    }

    /**
     * Gets a byte array representing the entire sequence. If the genome stores
     * differences rather than a complete sequence, then this may be an inefficient
     * way of accessing each state.
     *
     * @return
     */
	public byte[] getSequence() {
		byte[] sequence = new byte[masterSequence.length];
		System.arraycopy(masterSequence, 0, sequence, 0, masterSequence.length);
		for (int i = 0; i < mutations.length; i++) {
			sequence[mutations[i].position] = mutations[i].state;
		}
		return sequence;
	}

    public double getLogFitness() {
        return logFitness;
    }

    public void setLogFitness(double logFitness) {
        this.logFitness = logFitness;
    }

    public int getFrequency() {
        return frequency;
    }

    public void setFrequency(int frequency) {
        this.frequency = frequency;
    }

    // private members

    private double logFitness = 0.0;
    private int frequency = 1;

    private int mutationCount = 0;
    private Mutation[] mutations = null;

    // static members

    public static void setMasterSequence(byte[] masterSequence, double logFitness) {
        Genome.masterSequence = masterSequence;
        Genome.masterLogFitness = logFitness;
    }

    public static byte[] getMasterSequence() {
        return masterSequence;
    }

    public static double getMasterLogFitness() {
        return masterLogFitness;
    }

    private static byte[] masterSequence = null;
    private static double masterLogFitness = 0.0;
}