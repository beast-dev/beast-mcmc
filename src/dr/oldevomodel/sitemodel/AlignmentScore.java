/*
 * AlignmentScore.java
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

package dr.oldevomodel.sitemodel;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ExtractPairs;
import dr.evolution.alignment.GapUtils;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.math.DifferentialEvolution;
import dr.math.MultivariateFunction;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;

import java.io.FileReader;
import java.util.*;

/**
 * @author alexei
 * @version $Id: AlignmentScore.java,v 1.2 2005/04/11 11:29:37 alexei Exp $
 */
public class AlignmentScore implements UnivariateFunction, MultivariateFunction {

    SiteModel siteModel;
    ScoreMatrix scoreMatrix;
    SitePatterns sitePatterns;

    public AlignmentScore(ScoreMatrix scoreMatrix, SitePatterns sitePatterns) {
        this.scoreMatrix = scoreMatrix;
        this.siteModel = scoreMatrix.siteModel;
        this.sitePatterns = sitePatterns;
    }

    /**
     * compute function value
     *
     * @param argument function argument (vector)
     * @return function value
     */
    public double evaluate(double[] argument) {

        double kappa = argument[0];
        double time = argument[1];

        ((HKY)((GammaSiteModel)siteModel).getSubstitutionModel()).setKappa(kappa);

        scoreMatrix.setTime(time);

        // must return negative logLikelihood because optimizer finds minimum
        return -scoreMatrix.getScore(sitePatterns);
    }

    /**
     * compute function value
     *
     * @param time
     * @return function value
     */
    public double evaluate(double time) {
        scoreMatrix.setTime(time);

        // must return negative logLikelihood because optimizer finds minimum
        return -scoreMatrix.getScore(sitePatterns);
    }

    /**
     * get number of arguments
     *
     * @return number of arguments
     */
    public int getNumArguments() {
        return 2;
    }

    /**
     * get lower bound of argument
     *
     * @return lower bound
     */
    public double getLowerBound() {
        return 0.0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * get upper bound of argument
     *
     * @return upper bound
     */
    public double getUpperBound() {
        return 10.0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * get lower bound of argument n
     *
     * @param n argument number
     * @return lower bound
     */
    public double getLowerBound(int n) {

        return 0;
    }

    /**
     * get upper bound of argument n
     *
     * @param n argument number
     * @return upper bound
     */
    public double getUpperBound(int n) {
        if (n == 0) return 100;
        if (n == 1) return 100;
        return 10;
    }

    public static double[] getAlignmentScore(ScoreMatrix scoreMatrix, SitePatterns sitePatterns) {

        AlignmentScore alignmentScore = new AlignmentScore(scoreMatrix, sitePatterns);

        double[] params = new double[alignmentScore.getNumArguments()];
        DifferentialEvolution de = new DifferentialEvolution(params.length,params.length*10);
        for (int i = 0; i < params.length; i++) params[i] = 0.5;

        de.optimize(alignmentScore, params,1e-6, 1e-6);
        double score = alignmentScore.evaluate(params);

        double[] results = new double[params.length + 1];
        System.arraycopy(params, 0, results, 0, params.length);
        results[params.length] = score;

        return results;
    }

    public static double getGeneticDistance(ScoreMatrix scoreMatrix, SitePatterns sitePatterns) {

        AlignmentScore alignmentScore = new AlignmentScore(scoreMatrix, sitePatterns);

        double[] params = new double[alignmentScore.getNumArguments()];
        DifferentialEvolution de = new DifferentialEvolution(params.length,params.length*10);
        for (int i = 0; i < params.length; i++) params[i] = 0.5;

        de.optimize(alignmentScore, params,1e-6, 1e-6);
        //double score = alignmentScore.evaluate(params);

        return params[params.length-1];

    }

    public static double[] getFastAlignmentScore(ScoreMatrix scoreMatrix, SitePatterns sitePatterns) {

            AlignmentScore alignmentScore = new AlignmentScore(scoreMatrix, sitePatterns);

            UnivariateMinimum de = new UnivariateMinimum();

            double time = de.optimize(alignmentScore, 1e-6);
            double score = alignmentScore.evaluate(time);

        return new double[] {time, score};
        }


    private static void printFrequencyTable(List<Integer> list) {
        int max = 0;
        int total = 0;
        for (Integer size : list) {
            if (size > max) max = size;
            total += size;
        }
        int[] freqs = new int[max+1];
        for (Integer size : list) {
            freqs[size] += 1;
        }

        for (int i = 0; i < freqs.length; i++) {
            System.out.println(i + "\t" + freqs[i]);
        }
        System.out.println("Total = " + total);
    }

    public static void main(String[] args) throws java.io.IOException, Importer.ImportException {

        NexusImporter importer = new NexusImporter(new FileReader(args[0]));

        Alignment alignment = importer.importAlignment();

        ExtractPairs pairs = new ExtractPairs(alignment);

        Parameter muParam = new Parameter.Default(1.0);
        Parameter kappaParam = new Parameter.Default(1.0);
        kappaParam.addBounds(new Parameter.DefaultBounds(100.0, 0.0, 1));
        muParam.addBounds(new Parameter.DefaultBounds(1.0, 1.0, 1));

        Parameter freqParam = new Parameter.Default(alignment.getStateFrequencies());

        FrequencyModel freqModel = new FrequencyModel(Nucleotides.INSTANCE,freqParam);

        SubstitutionModel substModel = new HKY(kappaParam, freqModel);

        SiteModel siteModel = new GammaSiteModel(substModel, muParam, null, 1, null);

        ScoreMatrix scoreMatrix = new ScoreMatrix(siteModel, 0.1);


        double threshold = 0.1;

        List<PairDistance> pairDistances = new ArrayList<PairDistance>();
        Set<Integer> sequencesUsed = new HashSet<Integer>();
        List<Integer> allGaps = new ArrayList<Integer>();
        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            for (int j = i+1; j < alignment.getSequenceCount(); j++) {
                Alignment pairAlignment = pairs.getPairAlignment(i,j);
                if (pairAlignment != null) {
                    SitePatterns patterns = new SitePatterns(pairAlignment);

                    double distance = getGeneticDistance(scoreMatrix, patterns);
                    if (distance < threshold) {
                        List gaps = new ArrayList();
                        GapUtils.getGapSizes(pairAlignment, gaps);
                        pairDistances.add(new PairDistance(i,j,distance, gaps, pairAlignment.getSiteCount()));
                        System.out.print(".");
                    } else {
                        System.out.print("*");
                    }
                } else {
                    System.out.print("x");
                }
            }
            System.out.println();
        }

        Collections.sort(pairDistances);
        int totalLength = 0;
        for (PairDistance pairDistance : pairDistances) {

            Integer x = pairDistance.x;
            Integer y = pairDistance.y;

            if (!sequencesUsed.contains(x) && !sequencesUsed.contains(y)) {
                allGaps.addAll(pairDistance.gaps);
                sequencesUsed.add(x);
                sequencesUsed.add(y);
                System.out.println("Added pair (" + x + "," + y + ") d=" + pairDistance.distance + " L=" + pairDistance.alignmentLength);
                totalLength += pairDistance.alignmentLength;
            }
        }

        printFrequencyTable(allGaps);
        System.out.println("total length=" + totalLength);
    }
}
