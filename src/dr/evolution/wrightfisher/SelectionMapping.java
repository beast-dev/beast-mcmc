/*
 * SelectionMapping.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.wrightfisher;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

/**
 * @author Alexei Drummond
 *
 * @version $Id: SelectionMapping.java,v 1.2 2005/04/20 21:26:19 rambaut Exp $
 */
public class SelectionMapping {

   public static void main(String[] args) throws java.io.IOException {

		if (args.length != 1) throw new RuntimeException("Expect a file name containing parameters!");

		FileInputStream fi = new FileInputStream(args[0]);
		Properties props = new Properties();
		props.load(fi);
		fi.close();

		int genomeLength = Integer.parseInt(props.getProperty("genomeLength", "1000"));
		int populationSize = Integer.parseInt(props.getProperty("populationSize", "250"));
		int replicates = Integer.parseInt(props.getProperty("replicates", "1000"));
		double mu = Double.parseDouble(props.getProperty("mu", "1e-4"));
		double alpha = Double.parseDouble(props.getProperty("alpha", "0.06"));
        double pInv = Double.parseDouble(props.getProperty("pinv", "0.5"));
		double beta = Double.parseDouble(props.getProperty("beta", "1"));
		int stateSize = Integer.parseInt(props.getProperty("stateSize", "4"));
		int n = Integer.parseInt(props.getProperty("sampleSize", "40"));
		boolean randomFittest = Boolean.valueOf(props.getProperty("randomizeFittest", "false"));
		boolean outputAlignments = Boolean.valueOf(props.getProperty("outputAlignments", "false"));
		String alignmentFileName = props.getProperty("alignment.filename", "alignment.fasta");
		String fileName = props.getProperty("output.filename", "out.txt");
		int burninFactor = Integer.parseInt(props.getProperty("burninFactor", "10"));

		PrintWriter writer = new PrintWriter(new FileWriter(fileName));

		writer.println("// WMD v1.0");
		writer.println("// genomeLength: " + genomeLength);
		writer.println("// populationSize: " + populationSize);
		writer.println("// alpha: " + alpha);
		writer.println("// beta: " + beta);
		writer.println("// replicates: " + replicates);
		writer.println("// output.filename: " + fileName);
		writer.println("// mu: " + mu);
		writer.println("// outputAlignments: " + outputAlignments);
		writer.println("// randomFittest: " + randomFittest);
		writer.println("// alignment.filename: " + alignmentFileName);

		int[] totalMutationalDensity = new int[200];
		int[] totalLineages = new int[200];
		int[] nMutants = new int[genomeLength*3];
		double[] mrcaFitness = new double[1];

		final int generations = populationSize;

		ArrayList[] unfoldedSites = new ArrayList[populationSize+1];
		for (int i = 0; i < populationSize+1; i++) {
			unfoldedSites[i] = new ArrayList();
		}

		long start = System.currentTimeMillis();

		int dnaSingles = 0;
		int aaSingles = 0;
		int dnaSegs = 0;
		int aaSegs = 0;

		for (int reps = 0; reps < replicates; reps++) {
//			FitnessFunction f = null;
//			if (alpha == 0.0) {
//				f = new NeutralModel();
//			}
			FitnessFunction f = new CodonFitnessFunction(genomeLength, alpha, beta, pInv);
			Population p = new Population(populationSize, genomeLength*3, new SimpleMutator(mu, stateSize), f, randomFittest);

			Genome master = new SimpleGenome(genomeLength*3, f, randomFittest);

			p = Population.forwardSimulation(p, burninFactor * populationSize);

			int age = -1;
			while (age == -1) {
				p = Population.forwardSimulation(p, generations);
				age = p.getAgeOfMRCA(mrcaFitness);
			}
			Population children = Population.forwardSimulation(p, 1);
			writer.println(reps + "\t"+age + "\t" + p.getMeanParentFitness() + "\t" + mrcaFitness[0] + "\t" + p.getProportionAsFit(mrcaFitness[0]));

			if (reps % 10 == 0) {
				System.out.print(reps + "\t"+ age + "\t" + p.getMeanParentFitness() + "\t" + mrcaFitness[0] + "\t" + p.getProportionAsFit(mrcaFitness[0]));
				long millis = System.currentTimeMillis() - start;
				double seconds = millis / 1000.0;
				if (reps != 0) {
					double seconds2go = (replicates - reps) * seconds / reps;
					System.out.println(" -- "+ (Math.round(seconds2go / 36) / 100.0) + " hours");
				} else {
					System.out.println();
				}
			}

			//p.unfoldedSiteFrequencies(unfoldedSites);


			// ************************************************************************************
			// get the mutational density per lineage as a function of time
			// ************************************************************************************

			List<Integer>[] mutations = new ArrayList[200];
			for (int i = 0; i < mutations.length; i++) {
				mutations[i] = new ArrayList<Integer>(200);
			}

			int sampleSize = Math.min(n, populationSize);

			p.getMutationDensity(sampleSize, mutations);

			for (int i = 0; i < mutations.length; i++) {
				for (int j = 0; j < mutations[i].size(); j++) {
					totalMutationalDensity[i] += mutations[i].get(j);
				}
				totalLineages[i] += mutations[i].size();
			}

			// ************************************************************************************
			// output alignments if requested
			// ************************************************************************************

			SimpleAlignment alignment = new SimpleAlignment();
			SimpleAlignment aaAlignment = new SimpleAlignment();

			alignment.setDataType(Nucleotides.INSTANCE);
			aaAlignment.setDataType(AminoAcids.INSTANCE);

			if (outputAlignments) {
				PrintWriter pw = new PrintWriter(new FileWriter(alignmentFileName));


				//pw.println(sampleSize + "\t" + p.getGenomeLength());
				for (int i = 0; i < sampleSize; i++) {
					pw.print(">sequence_");
					if (i < 10) {
						pw.print("0");
					}
					pw.println(""+i);
					String dnaSequence = p.getGenome(i).getDNASequenceString();
					String aaSequence = p.getGenome(i).getAminoAcidSequenceString();
					pw.println(dnaSequence);
					alignment.addSequence(new Sequence(dnaSequence));
					aaAlignment.addSequence(new Sequence(aaSequence));
				}
				pw.println();
				pw.close();
			}

			// ************************************************************************************
			// calculate hamming distances from master sequence
			// ************************************************************************************
			for (int i = 0; i < sampleSize; i++) {
				nMutants[master.hammingDistance(p.getGenome(i))] += 1;
			}

			dnaSingles += countSingletons(alignment);
			aaSingles += countSingletons(aaAlignment);
			dnaSegs += countSegregating(alignment);
			aaSegs += countSegregating(aaAlignment);
		}




		for (int i = 0; i < unfoldedSites.length; i++) {
			double meanFitness = 0.0;
			double proportionNeutral = 0.0;
			double proportionPositive = 0.0;
			double proportionNegative = 0.0;
			for (int j = 0; j < unfoldedSites[i].size(); j++) {
				double relativeFitness = (Double) unfoldedSites[i].get(j);
				meanFitness += relativeFitness;
				if (relativeFitness == 1.0) {
					proportionNeutral += 1.0;
				} else if (relativeFitness > 1.0) {
					proportionPositive += 1.0;
				} else {
					proportionNegative += 1.0;
				}

			}
			meanFitness /= (double)unfoldedSites[i].size();
			proportionNeutral /= (double)unfoldedSites[i].size();
			proportionPositive /= (double)unfoldedSites[i].size();
			proportionNegative /= (double)unfoldedSites[i].size();
			writer.println(i + "\t" + unfoldedSites[i].size() +
				"\t" + meanFitness + "\t" + proportionNeutral + "\t" + proportionPositive + "\t" + proportionNegative);

		}


		writer.println("--------------------");
		writer.println("SINGLETON COUNTS");
		writer.println("--------------------");
		writer.println("dna singletons = " + dnaSingles);
		writer.println("aa singletons = " + aaSingles);

		int dnaNons = dnaSegs-dnaSingles;
		int aaNons = aaSegs-aaSingles;

		System.out.println("dna singletons = " + dnaSingles);
		System.out.println("aa singletons = " + aaSingles);
		System.out.println("dna segregating = " + dnaSegs);
		System.out.println("aa segregating = " + aaSegs);
		System.out.println("dna non-singles = " + dnaNons);
		System.out.println("aa non-singles = " + aaNons);
		System.out.println("ratio = " + ((double)dnaSingles/(double)aaSingles));
		System.out.println("ratio(non) = " + ((double)dnaNons/(double)aaNons));

		writer.println("--------------------");
		writer.println("MUTATIONAL DENSITIES");
		writer.println("--------------------");
		for (int i = 0; i < 200; i++) {
			writer.println(totalMutationalDensity[i] + "\t" + totalLineages[i]);
		}

		// ************************************************************************************
		// output hamming distances
		// ************************************************************************************
		writer.println("--------------------");
		writer.println("Hamming distance distribution");
		writer.println("--------------------");
		for (int i = 0; i < nMutants.length; i++) {
			writer.println(i + "\t" + nMutants[i]);
		}
		writer.close();
	}

	public static int countSingletons(Alignment a) {

		int count = 0;
		int[] counts = new int[a.getDataType().getStateCount()];
		for (int i = 0; i < a.getSiteCount(); i++) {
			count += (isSingleton(a, i, counts)?1:0);
		}

		return count;
	}

	public static int countSegregating(Alignment a) {

		int count = 0;
		int[] counts = new int[a.getDataType().getStateCount()];
		for (int i = 0; i < a.getSiteCount(); i++) {
			count += (isSegregating(a, i, counts)?1:0);
		}

		return count;
	}

	public static boolean isSingleton(Alignment a, int j, int[] counts) {
		for (int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
		for (int i = 0; i < a.getSequenceCount(); i++) {
			int state = a.getState(i, j);
			if (state >= 0 && state < counts.length) {
				counts[a.getState(i, j)] += 1;
			}
		}
        for(int count : counts) {
            if( count == 1 ) return true;
        }
		return false;
	}

	public static boolean isSegregating(Alignment a, int j, int[] counts) {
		for (int i = 0; i < counts.length; i++) {
			counts[i] = 0;
		}
		int seqCount = a.getSequenceCount();
		for (int i = 0; i < seqCount; i++) {
			int state = a.getState(i, j);
			if (state >= 0 && state < counts.length) {
				counts[a.getState(i, j)] += 1;
			}
		}
        for(int count : counts) {
            if( count > 0 && count < seqCount ) return true;
        }
		return false;
	}

}
