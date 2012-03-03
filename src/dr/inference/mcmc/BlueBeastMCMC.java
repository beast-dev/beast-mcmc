//package dr.inference.mcmc;
//
//import bb.main.BlueBeast;
//import dr.inference.loggers.Logger;
//import dr.inference.model.Likelihood;
//import dr.inference.operators.MCMCOperator;
//import dr.inference.operators.OperatorSchedule;
//import dr.inference.operators.SimpleOperatorSchedule;
//import dr.inference.prior.Prior;
//
//import bb.main.BlueBeast;
//import dr.inference.markovchain.BlueBeastMarkovChain;
//import dr.inference.loggers.Logger;
//import dr.inference.mcmc.MCMCCriterion;
//import dr.inference.mcmc.MCMCOptions;
//import dr.inference.model.Likelihood;
//import dr.inference.operators.*;
//import dr.inference.operators.OperatorSchedule;
//import dr.inference.prior.Prior;
//
///**
// * Created by IntelliJ IDEA.
// * User: sibon
// * Date: 3/15/11
// * Time: 2:43 PM
// * To change this template use File | Settings | File Templates.
// *
// * @author Wai Lok Sibon Li
// *
// */
//
////@Description("Sibon's MCMC extension with periodical tests")
//public class BlueBeastMCMC extends MCMC {
//
//    public BlueBeastMCMC(String s) {
//        super(s);
//    }
//
//
//    /**
//     * Must be called before calling chain.
//     *
//     * @param options    the options for this MCMC analysis
//     * @param schedule   operator schedule to be used in chain.
//     * @param likelihood the likelihood for this MCMC
//     * @param loggers    an array of loggers to record output of this MCMC run
//     * @param bb         a Blue Beast control object
//     */
//    public void init(
//            MCMCOptions options,
//            Likelihood likelihood,
//            OperatorSchedule schedule,
//            Logger[] loggers, BlueBeast bb) {
//
//        init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers, bb);
//    }
//
//    /**
//     * Must be called before calling chain.
//     *
//     * Override method from MCMC.java
//     *
//     * @param options    the options for this MCMC analysis
//     * @param prior      the prior disitrbution on the model parameters.
//     * @param schedule   operator schedule to be used in chain.
//     * @param likelihood the likelihood for this MCMC
//     * @param loggers    an array of loggers to record output of this MCMC run
//     * @param bb         a Blue Beast control object
//     */
//    public void init(
//            MCMCOptions options,
//            Likelihood likelihood,
//            Prior prior,
//            OperatorSchedule schedule,
//            Logger[] loggers, BlueBeast bb) {
//
//
//        MCMCCriterion criterion = new MCMCCriterion();
//        criterion.setTemperature(options.getTemperature());
//
////        mc = new BlueBeastMarkovChain(prior, likelihood, schedule, criterion,
////                options.fullEvaluationCount(), options.minOperatorCountForFullEvaluation(), options.useCoercion());
//        mc = new BlueBeastMarkovChain(prior, likelihood, schedule, criterion,
//                options.fullEvaluationCount(), options.minOperatorCountForFullEvaluation(), options.useCoercion(), bb);
//
//        this.options = options;
//        this.loggers = loggers;
//        this.schedule = schedule;
//        //this.init(options, likelihood, prior, schedule, loggers);
//
//        //initialize transients
//        currentState = 0;
//        this.bb = bb;
//    }
//
//
//    /**
//     * Must be called before calling chain.
//     *
//     * @param chainLength chain length
//     * @param likelihood the likelihood for this MCMC
//     * @param operators  an array of MCMC operators
//     * @param loggers    an array of loggers to record output of this MCMC run
//     * @param bb         a Blue Beast control object
//     */
//    public void init(long chainLength,
//                     Likelihood likelihood,
//                     MCMCOperator[] operators,
//                     Logger[] loggers, BlueBeast bb) {
//
//        MCMCOptions options = new MCMCOptions();
//        options.setCoercionDelay(0);
//        options.setChainLength(chainLength);
//        MCMCCriterion criterion = new MCMCCriterion();
//        criterion.setTemperature(1);
//        OperatorSchedule schedule = new SimpleOperatorSchedule();
//        for (MCMCOperator operator : operators) schedule.addOperator(operator);
//
//        init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers, bb);
//
//        this.bb = bb;
//    }
//
//
//    public long getCheckInterval() {
//        return checkInterval;
//    }
//
//    public void setCheckInterval(long checkInterval) {
//        this.checkInterval = checkInterval;
//    }
//
//
//
//    /**
//     * this markov chain does most of the work.
//     */
//    //private MarkovChain mc;
////    private BlueBeastMarkovChain mc;
////    private MCMCOptions options;                   // merge
////    private Logger[] loggers;
////    private OperatorSchedule schedule;
////    private String id = null;
//
//    private long checkInterval; /* period between samples being tested (default 1000) */
//    //private boolean doCheck; /* flag to indicate whether samples should be tested (default true) */
//
//
//
//
//    /**
//     * This method actually initiates the MCMC analysis.
//     * This is where the BlueBeastMCMC magic happens
//     */
//    public void chain() {
//
//        stopping = false;
//        currentState = 0;
//
//        timer.start();
//
//        if (loggers != null) {
//            for (Logger logger : loggers) {
//                logger.startLogging();
//            }
//        }
//
//        if (!stopping) {
//            mc.addMarkovChainListener(chainListener);
//
//
//            long chainLength = getChainLength();
//
//            final int coercionDelay = getCoercionDelay();
//
//            if (coercionDelay > 0) {
//                // Run the chain for coercionDelay steps with coercion disabled
//                mc.runChain(coercionDelay, true);
//                chainLength -= coercionDelay;
//
//                // reset operator acceptance levels
//                for (int i = 0; i < schedule.getOperatorCount(); i++) {
//                    schedule.getOperator(i).reset();
//                }
//            }
//
//            //mc.runChain(chainLength, false);
//            mc.runChain(bb.getMaxChainLength(), false);
//
//            mc.terminateChain();
//
//            mc.removeMarkovChainListener(chainListener);
//        }
//        timer.stop();
//    }
//
//    private BlueBeast bb;
//
//}
//
