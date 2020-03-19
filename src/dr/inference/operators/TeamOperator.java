/*
 * TeamOperator.java
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

package dr.inference.operators;

import dr.math.MathUtils;

/**
 * Combines N oprators - each time a random subset of k out of the N are run. If operators have weight then the
 * probability of any subset is proportional to the sum of it's members weights.
 *
 * Currently very basic as no Coercible operators are handled correctly. Getting coercion to work is non trivial.
 *
 * @author  Joseph Heled
 */

public class TeamOperator extends SimpleMCMCOperator /*implements AdaptableMCMCOperator*/ {

    private final MCMCOperator[] operators;
    //private final ArrayList<Integer> operatorToOptimizeList;

    //private int currentOptimizedOperator;
   // private final double targetProbability;
    private final int nPick;
    private final boolean unequalWeights;

    private final MCMCOperator[] currentRound;
    private int nToReject;
    private int[][] binomial;

    public TeamOperator(MCMCOperator[] operators, int nPick, double weight) {
        setWeight(weight);

        this.operators = operators;
        //targetProbability = targetProb;
        final int N = operators.length;
        assert 0 < nPick && nPick <= N;
        this.nPick = nPick;
        currentRound = new MCMCOperator[N];


        {
            boolean b = false;
            double w = operators[0].getWeight();
            for( MCMCOperator o : operators ) {
                if( o.getWeight() != w ) {
                    b = true;
                    break;
                }
            }
            unequalWeights = b;
        }

        if( unequalWeights ) {
            final int M = N+1;
            binomial = new int[M][M];
            for(int n = 0; n < M; ++n) {
        //        binomial[n] = new int[operators.length];
              binomial[n][0] = 1;
            }

            for(int n = 1; n < M; ++n) {
              for(int k = 1; k <= n; ++k) {
                    binomial[n][k] = binomial[n-1][k] + binomial[n-1][k-1];
                }
            }
        } else {
            for(int k = 0; k < N; ++k) {
                currentRound[k] = operators[k];
            }
        }
    }

//    public void addOperator(SimpleMCMCOperator operation) {
//
//        operatorList.add(operation);
//        if (operation instanceof AdaptableMCMCOperator) {
//
//            if (((AdaptableMCMCOperator) operation).getMode() == AdaptationMode.ADAPTATION_ON)
//
//                operatorToOptimizeList.add(operatorList.size() - 1);
//
//        }
//    }

    private void choose() {
        final int n = operators.length;
        if( nPick < n ) {
            if( unequalWeights ) {
              chooseUsingWeights();
            } else {
                // equal weights, just pick a subset of 'nPick' operators uniformly
                for(int k = 0; k < nPick; ++k) {
                    final int which = k + MathUtils.nextInt(n - k);
                    final MCMCOperator tmp = currentRound[k];
                    currentRound[k] = currentRound[which];
                    currentRound[which] = tmp;
                }
            }
        }
    }

    private void chooseUsingWeights() {
        // sum(o_w : o in operators already selected)
        double inSumWeights = 0.0;
        // sum(o_w : o in remaining operators)
        double sumWeightsRemaining = 0;
        for( MCMCOperator o : operators ) {
            sumWeightsRemaining += o.getWeight();
        }

        // Number of operators still to pick
        int k = nPick;
        // Operator under consideration
        int j = 0;
        while( k > 0 ) {
            // remaining to choose from
            final int n = operators.length - j;

            if( k == n ) {
                // speedup
                for( ; k > 0; k--, j++ ) {
                    currentRound[k-1] = operators[j];
                }
            } else {
                final int cnk = binomial[n][k];
                final int cnk1 = binomial[n-1][k-1];
                final int cnk2 = k >= 2 ? binomial[n-2][k-2] : 0;

                final double tot = cnk1 * sumWeightsRemaining +  cnk * inSumWeights;
                final double we0 = operators[j].getWeight();
                final double has = cnk2 * (sumWeightsRemaining-we0) + cnk1 * (we0+inSumWeights);

                final double r = MathUtils.nextDouble();
                if( r < has/tot ) {
                    currentRound[k-1] = operators[j];
                    k -= 1;
                    inSumWeights += we0;
                }
                j += 1;
                sumWeightsRemaining -= we0;
            }
         }
     }

    public final double doOperation() {

        choose();

        double logP = 0;

        for(int k = 0; k < nPick; ++k) {
            MCMCOperator operation = currentRound[k];
            logP += operation.operate();
        }
        nToReject = nPick;

        return logP;
    }

    public void accept(double deviation) {
        super.accept(deviation);
        for(int k = 0; k < nPick; ++k) {
           currentRound[k].accept(deviation);
        }
    }

    public void reject() {
       super.reject();
        for(int k = 0; k < nToReject; ++k) {
           currentRound[k].reject();
        }
    }

    public void reset() {
        for( MCMCOperator op : operators ) {
            op.reset();
        }
    }

    public String getOperatorName() {
       StringBuffer sb = new StringBuffer("Team " + nPick + " (");
       for( MCMCOperator operation : operators ) {
           sb.append(operation.getOperatorName()+",");
       }
       return sb.substring(0, sb.length()-1) + ")";
    }

//    public double getTargetAcceptanceProbability() {
//        return targetProbability;
//    }

    public String getPerformanceSuggestion() {
        return "";
    }

}