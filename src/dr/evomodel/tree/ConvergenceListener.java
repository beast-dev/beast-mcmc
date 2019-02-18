/*
 * ConvergenceListener.java
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

/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.tree.Tree;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Model;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * @author shhn001
 */
@Deprecated
public class ConvergenceListener implements MarkovChainListener {

    private MarkovChain chain;

    private final int LOG_EVERY = 10;

    private String outputFilename;

    private Convergence convergence;

    private double threshold;

    private long startTime;

    private double distance;

    /**
     *
     */
    public ConvergenceListener(MarkovChain chain, String cladeFilename, String outputfile, Tree tree, double threshold) {
        convergence = new Convergence(tree, LOG_EVERY, cladeFilename);

        this.chain = chain;

        this.outputFilename = outputfile;
        File f = new File(outputFilename);
        if (f.exists()) {
            f.delete();
        }

        this.threshold = threshold;

        startTime = System.currentTimeMillis();
    }

    /* (non-Javadoc)
      * @see dr.inference.markovchain.MarkovChainListener#bestState(int, dr.inference.model.Model)
      */
    public void bestState(long state, MarkovChain markovChain, Model bestModel) {
        // do nothing
    }

    /* (non-Javadoc)
      * @see dr.inference.markovchain.MarkovChainListener#currentState(int, dr.inference.model.Model)
      */
    public void currentState(long state, MarkovChain markovChain, Model currentModel) {
        distance = convergence.log(state);

        if (distance <= threshold) {
            chain.pleaseStop();
        }
    }

    /* (non-Javadoc)
      * @see dr.inference.markovchain.MarkovChainListener#finished(int)
      */
    public void finished(long chainLength, MarkovChain markovChain) {
        // write the time used
        long time = System.currentTimeMillis() - startTime;

        try {
            File f = new File(outputFilename);
            FileWriter fw = new FileWriter(f);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("" + time);
            bw.newLine();
            bw.write("" + distance);
            bw.newLine();
            bw.write("" + chainLength);
            bw.newLine();
            bw.close();
            fw.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        System.out.println("Time needed:\t\t" + time);
        System.out.println("Iterations needed:\t" + chainLength);
        System.out.println("Time per Iteration:\t" + ((double) time) / chainLength);
    }

}
