package dr.inference.checkpoint;

import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainDelegate;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.operators.OperatorSchedule;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CheckPointMarkovChainDelegate implements MarkovChainDelegate {
    public final static String LOAD_DUMP_FILE = "load.dump.file";
    public final static String SAVE_DUMP_FILE = "save.dump.file";
    public final static String DUMP_STATE = "dump.state";
    public final static String DUMP_EVERY = "dump.every";

    private String dumpStateFile = null;

    @Override
    public void setup(MCMCOptions options, OperatorSchedule schedule, MarkovChain markovChain) {
        dumpStateFile = System.getProperty(LOAD_DUMP_FILE);
        String fileName = System.getProperty(SAVE_DUMP_FILE, null);
        if (System.getProperty(DUMP_STATE) != null) {
            long debugWriteState = Long.parseLong(System.getProperty(DUMP_STATE));
            markovChain.addMarkovChainListener(new CheckPointChainListener(markovChain, debugWriteState, false, fileName));
        }
        if (System.getProperty(DUMP_EVERY) != null) {
            long debugWriteEvery = Long.parseLong(System.getProperty(DUMP_EVERY));
            markovChain.addMarkovChainListener(new CheckPointChainListener(markovChain, debugWriteEvery, true, fileName));
        }
    }

    @Override
    public void currentState(long state) {

    }

    @Override
    public void currentStateEnd(long state) {

    }

    @Override
    public void finished(long chainLength) {

    }
}
