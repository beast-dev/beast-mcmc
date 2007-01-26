package dr.inference.parallel;

import dr.inference.loggers.Logger;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.AbstractModel;
import dr.util.Identifiable;
import dr.xml.*;
import mpi.MPI;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 23, 2007
 * Time: 9:47:41 AM
 * To change this template use File | Settings | File Templates.
 */
public class MPILikelihoodRunner implements Runnable, Identifiable {

    public static final String PARALLEL_CALCULATOR = "parallelCalculator";
    public static final int origin = 0;

    public MPILikelihoodRunner(Likelihood likelihood) {
       this.likelihood = likelihood;
       this.model = likelihood.getModel();
    }

    public void run() {
        System.err.println("Calculator is running");
        System.err.println("rank = " + mpiRank);
        System.err.println("size = " + mpiSize);

        while (!terminate) {

            // Wait for calculation request
            int[] jobId = new int[1];
             MPI.COMM_WORLD.Recv(jobId,0,1,MPI.INT,origin,ServiceRequest.MSG_REQUEST_TYPE);
            

        }



        finalize();
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;

    }

    /**
     * Initialize connection to MPI world
     */
    public void init() {
        //  MPI.Init(null);
        mpiRank = MPI.COMM_WORLD.Rank();
        mpiSize = MPI.COMM_WORLD.Size();
    }

    public void finalize() {
        MPI.Finalize();
    }

    /**
     * XML Parser
     */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARALLEL_CALCULATOR;
        }

        /**
         * @return a tree object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

           // runner.init();
            Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
            ArrayList loggers = new ArrayList();

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Logger) {
                    loggers.add(child);
                }
            }

            Logger[] loggerArray = new Logger[loggers.size()];
            loggers.toArray(loggerArray);

            java.util.logging.Logger.getLogger("dr.inference").info("Creating the parallelCalculator chain:");
            /*+
        "\n  chainLength=" + options.getChainLength() +
        "\n  autoOptimize=" + options.useCoercion());*/

            MPILikelihoodRunner runner = new MPILikelihoodRunner(likelihood);

            runner.init();

            return runner;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an MPI-based likelihood calculator.";
        }

        public Class getReturnType() {
            return MPILikelihoodRunner.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//                AttributeRule.newIntegerRule(CHAIN_LENGTH),
//                AttributeRule.newBooleanRule(COERCION, true),
//                AttributeRule.newIntegerRule(PRE_BURNIN, true),
//                new ElementRule(OperatorSchedule.class ),
                new ElementRule(Likelihood.class),
//                new ElementRule(Logger.class, 1, Integer.MAX_VALUE )
        };

    };

    // Private variables

    private String id;
    private int mpiRank;
    private int mpiSize;
    private boolean terminate = false;
    private Likelihood likelihood;
    private Model model;

}
