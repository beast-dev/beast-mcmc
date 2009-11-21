package dr.app.beagle.evomodel.substmodel;

import dr.inference.markovjumps.MarkovJumpsCore;
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
        super(substModel.getModelName());
        this.substModel = substModel;
        this.eigenDecomposition = substModel.getEigenDecomposition();
        stateCount = substModel.getDataType().getStateCount();
        markovJumpsCore = new MarkovJumpsCore(stateCount);
        setupStorage();
        addModel(substModel);
    }

    protected void setupStorage() {
        rateMatrix = new double[stateCount * stateCount];
        transitionProbs = new double[stateCount * stateCount];
        rateReg = new double[stateCount * stateCount];
    }

    private void makeRateRegistrationMatrix(double[] registration,
                                            double[] rateReg) {

        substModel.getInfinitesimalMatrix(rateMatrix);
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                rateReg[index] = rateMatrix[index] * registration[index];
                index++;
            }
        }
    }

    public void computeCondMeanMarkovJumps(double[] registration,
                                              double time,
                                              double[] countMatrix) {

        substModel.getTransitionProbabilities(time, transitionProbs);

        // TODO this requires work, only do once per traversal
        makeRateRegistrationMatrix(registration, rateReg);

        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        markovJumpsCore.computeCondMeanMarkovJumps(evec, ievc, eval, rateReg, time, transitionProbs, countMatrix);
    }

    public void computeJointMeanMarkovJumps(double[] registration,
                                               double time,
                                               double[] countMatrix) {

        // TODO this requires work, only do once per traversal
        makeRateRegistrationMatrix(registration, rateReg);

        double[] evec = eigenDecomposition.getEigenVectors();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();

        markovJumpsCore.computeJointMeanMarkovJumps(evec, ievc, eval, rateReg, time, countMatrix);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
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
        
        MarkovJumpsSubstitutionModel markovjumps = new MarkovJumpsSubstitutionModel(substModel);
        double[] r = new double[states * states];
        double[] q = new double[states * states];
        double[] j = new double[states * states];
        double[] c = new double[states * states];
        double[] p = new double[states * states];

        double time = 1.0;
        int from = 0; // A
        int to = 1; // C
        MarkovJumpsCore.fillRegistrationMatrix(r, from, to, states);

        substModel.getInfinitesimalMatrix(q);

        substModel.getTransitionProbabilities(time, p);

        markovjumps.computeJointMeanMarkovJumps(r, time, j);

        markovjumps.computeCondMeanMarkovJumps(r, time, c);

        MarkovJumpsCore.makeComparableToRPackage(q);
        System.err.println("Q = " + new Vector(q));

        MarkovJumpsCore.makeComparableToRPackage(r);
        System.err.println("R = " + new Vector(r));

        MarkovJumpsCore.makeComparableToRPackage(p);
        System.err.println("P = " + new Vector(p));

        MarkovJumpsCore.makeComparableToRPackage(j);
        System.err.println("J = " + new Vector(j));

        MarkovJumpsCore.makeComparableToRPackage(c);
        System.err.println("C = " + new Vector(c));
    }

    private int stateCount;
    private double[] rateReg;
    private double[] transitionProbs;
    private double[] rateMatrix;

    private SubstitutionModel substModel;
    private EigenDecomposition eigenDecomposition;
    private MarkovJumpsCore markovJumpsCore;
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
