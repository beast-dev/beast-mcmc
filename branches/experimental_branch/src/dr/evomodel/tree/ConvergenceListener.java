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
    public void bestState(long state, Model bestModel) {
        // do nothing
    }

    /* (non-Javadoc)
      * @see dr.inference.markovchain.MarkovChainListener#currentState(int, dr.inference.model.Model)
      */
    public void currentState(long state, Model currentModel) {
        distance = convergence.log(state);

        if (distance <= threshold) {
            chain.pleaseStop();
        }
    }

    /* (non-Javadoc)
      * @see dr.inference.markovchain.MarkovChainListener#finished(int)
      */
    public void finished(long chainLength) {
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
