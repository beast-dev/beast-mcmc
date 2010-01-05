package dr.app.beagle.evomodel.substmodel;

import dr.inference.markovjumps.MarkovJumpsCore;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.evolution.datatype.Nucleotides;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A base class for implementing Markov chain-induced counting processes (markovjumps) in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 *         Journal of Mathematical Biology, 56, 391-412.
 */

public class MarkovJumpsSubstitutionModel extends AbstractModel {

    public MarkovJumpsSubstitutionModel(SubstitutionModel substModel) {
        this(substModel, MarkovJumpsType.COUNTS);
    }

    public MarkovJumpsSubstitutionModel(SubstitutionModel substModel, MarkovJumpsType type) {
        super(substModel.getModelName());
        this.substModel = substModel;
        this.eigenDecomposition = substModel.getEigenDecomposition();
        stateCount = substModel.getDataType().getStateCount();
        markovJumpsCore = new MarkovJumpsCore(stateCount);
        this.type = type;
        setupStorage();
        addModel(substModel);
    }

    protected void setupStorage() {
        rateMatrix = new double[stateCount * stateCount];
        transitionProbs = new double[stateCount * stateCount];
        rateReg = new double[stateCount * stateCount];
        registration = new double[stateCount * stateCount];
    }

    public void setRegistration(double[] inRegistration) {

        if (type == MarkovJumpsType.COUNTS) {
                   
            System.arraycopy(inRegistration, 0, registration, 0, stateCount * stateCount);
            for (int i = 0; i < stateCount; i++) {
                registration[i * stateCount + i] = 0;  // diagonals are zero
            }

        } else if (type == MarkovJumpsType.REWARDS) {

            int index = 0;
            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    if (i == j) {
                        registration[index] = inRegistration[i];
                    } else {
                        registration[index] = 0; // Off-diagonals are zero
                    }
                    index++;
                }
            }

        } else {
            throw new RuntimeException("Unknown expectation type in MarkovJumps");
        }
        regRateChanged = true;
    }

    private void makeRateRegistrationMatrix(double[] registration,
                                            double[] rateReg) {

        if (type == MarkovJumpsType.COUNTS) {

            substModel.getInfinitesimalMatrix(rateMatrix);
            int index = 0;
            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    rateReg[index] = rateMatrix[index] * registration[index];
                    index++;
                }
            }

        } else if (type == MarkovJumpsType.REWARDS) {

            System.arraycopy(registration,0,rateReg,0,stateCount * stateCount);

        } else {
            throw new RuntimeException("Unknown expectation type in MarkovJumps");
        }

        regRateChanged = false;
    }

    public void computeCondMeanMarkovJumps(double time,
                                           double[] countMatrix) {

        substModel.getTransitionProbabilities(time, transitionProbs);
        computeCondMeanMarkovJumps(time,transitionProbs,countMatrix);
    }

    public void computeCondMeanMarkovJumps(double time,
                                           double[] transitionProbs,
                                           double[] countMatrix) {

        if (regRateChanged) {
            makeRateRegistrationMatrix(registration,rateReg);
        }

        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        markovJumpsCore.computeCondMeanMarkovJumps(evec, ievc, eval, rateReg, time, transitionProbs, countMatrix);
    }

    public void computeJointMeanMarkovJumps(double time,
                                            double[] countMatrix) {

        if (regRateChanged) {
            makeRateRegistrationMatrix(registration,rateReg);
        }

        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        markovJumpsCore.computeJointMeanMarkovJumps(evec, ievc, eval, rateReg, time, countMatrix);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == substModel && type == MarkovJumpsType.COUNTS) {
            regRateChanged = true;
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    protected void storeState() {
        // Do nothing
    }

    protected void restoreState() {
        // Do nothing
    }

    protected void acceptState() {
        // Do nothing
    }


    public static void main(String[] args) {

        HKY substModel = new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.3, 0.2, 0.25, 0.25})); // A,C,G,T
        int states = substModel.getDataType().getStateCount();
        
        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(substModel,
                MarkovJumpsType.COUNTS);
        double[] r = new double[states * states];
        double[] q = new double[states * states];
        double[] j = new double[states * states];
        double[] c = new double[states * states];
        double[] p = new double[states * states];

        double time = 1.0;
        int from = 0; // A
        int to = 1; // C
        MarkovJumpsCore.fillRegistrationMatrix(r, from, to, states, 1.0);
        markovjumps.setRegistration(r);

        substModel.getInfinitesimalMatrix(q);

        substModel.getTransitionProbabilities(time, p);

        markovjumps.computeJointMeanMarkovJumps(time, j);

        markovjumps.computeCondMeanMarkovJumps(time, c);

        MarkovJumpsCore.makeComparableToRPackage(q);
        System.err.println("Q = " + new Vector(q));

        MarkovJumpsCore.makeComparableToRPackage(p);
        System.err.println("P = " + new Vector(p));

        System.err.println("Counts:");
        MarkovJumpsCore.makeComparableToRPackage(r);
        System.err.println("R = " + new Vector(r));

        MarkovJumpsCore.makeComparableToRPackage(j);
        System.err.println("J = " + new Vector(j));

        MarkovJumpsCore.makeComparableToRPackage(c);
        System.err.println("C = " + new Vector(c));

        markovjumps = new MarkovJumpsSubstitutionModel(substModel, MarkovJumpsType.REWARDS);
        double[] rewards = {1.0, 1.0, 1.0, 1.0};
        markovjumps.setRegistration(rewards);

        markovjumps.computeJointMeanMarkovJumps(time, j);

        markovjumps.computeCondMeanMarkovJumps(time, c);

        System.err.println("Rewards:");
        MarkovJumpsCore.makeComparableToRPackage(rewards);
        System.err.println("R = " + new Vector(rewards));

        MarkovJumpsCore.makeComparableToRPackage(j);
        System.err.println("J = " + new Vector(j));

        MarkovJumpsCore.makeComparableToRPackage(c);
        System.err.println("C = " + new Vector(c));                        
    }

    private int stateCount;
    private double[] rateReg;
    private double[] transitionProbs;
    private double[] rateMatrix;
    private double[] registration;

    private SubstitutionModel substModel;
    private EigenDecomposition eigenDecomposition;
    private MarkovJumpsCore markovJumpsCore;

    private boolean regRateChanged = true;

    private MarkovJumpsType type;
}

/*
# R script to compare main to original package:

subst.model.eigen = as.eigen.hky(c(2,1),
	c(0.3,0.25,0.2,0.25),scale=T)

time = 1.0
from = 0 # A
to   = 2 # C
R = matrix(0,nrow=4,ncol=4)
R[from + to*4 + 1] = 1.0

P = matexp.eigen(subst.model.eigen,time)
J = joint.mean.markov.jumps(subst.model.eigen,R,time)
C = cond.mean.markov.jumps(subst.model.eigen,R,time)

 */
