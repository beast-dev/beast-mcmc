package dr.inference.model;


import org.apache.commons.math.special.Beta;

/**
 * @author Max Tolkoff
 */
public class IndianBuffetProcessPrior extends AbstractModelLikelihood {

    public IndianBuffetProcessPrior(Parameter alpha, Parameter beta, MatrixParameter data) {
        super(null);
        this.alpha=alpha;
        alpha.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0, 1));
        addVariable(alpha);
        this.beta=beta;
        beta.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0, 1));
        addVariable(beta);
        this.data=data;
        addVariable(data);
    }

    private int factorial(int num){
        if(num<0){
            throw new RuntimeException("Cannot take a negative factorial");
        }
        else if(num==0){
            return 1;
        }
        else
        {
            int fac=1;
            for (int i = 0; i <num ; i++) {
                fac*=(i+1);
            }
            return fac;
        }
    }

    private double H(){
        if(!betaKnown) {
            H = 0;
            for (int i = 0; i < data.getRowDimension(); i++) {
                H += beta.getParameterValue(0) / (beta.getParameterValue(0) + i);
            }
        }
        return H;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown=false;
        if(variable==beta)
            betaKnown=false;
        if(variable==data)
            dataKnown=false;

    }

    @Override
    protected void storeState() {
        storedBetaKnown=betaKnown;
        storedContainsNonZeroElements=containsNonZeroElements;
        storedDataKnown=dataKnown;
        storedLikelihoodKnown=likelihoodKnown;
        storedLogLikelihood=logLikelihood;
        storedRowCount=rowCount;
        storedKPlus=KPlus;
        storedH=H;
        storedBottom=bottom;
        storedSum2=sum2;

    }

    @Override
    protected void restoreState() {
        betaKnown=storedBetaKnown;
        containsNonZeroElements=storedContainsNonZeroElements;
        dataKnown=storedDataKnown;
        likelihoodKnown=storedLikelihoodKnown;
        logLikelihood=storedLogLikelihood;
        rowCount=storedRowCount;
        KPlus=storedKPlus;
        H=storedH;
        bottom=storedBottom;
        sum2=storedSum2;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if(!likelihoodKnown){
        logLikelihood=calculateLogLikelihood();
            likelihoodKnown=true;
        }
        return logLikelihood;
    }

    private double calculateLogLikelihood(){

        int sum;

        if(!dataKnown) {
            bottom=1;
            boolean[] isExplored= new boolean[data.getColumnDimension()];
            containsNonZeroElements = new boolean[data.getColumnDimension()];
            rowCount = new int[data.getColumnDimension()];
            boolean same;
            for (int i = 0; i < data.getRowDimension(); i++) {
                if (data.getParameterValue(i, 0) != 0)
                    containsNonZeroElements[0] = true;
            }
            for (int i = 0; i < data.getColumnDimension(); i++) {
                sum = 1;
                if (!isExplored[i]) {
                    for (int j = i + 1; j < data.getColumnDimension(); j++) {
                        same = true;
                        if (!isExplored[j]) {
                            for (int k = 0; k < data.getRowDimension(); k++) {
                                if (data.getParameterValue(k, i) != data.getParameterValue(k, j))
                                    same = false;
                                if (data.getParameterValue(k, j) != 0) {
                                    containsNonZeroElements[j] = true;
                                }
//                        rowCount[j]+=data.getParameterValue(k,j);
                            }
                        }
                        if (same && containsNonZeroElements[j]) {
                            isExplored[j] = true;
                            sum += 1;
                        } else if (!containsNonZeroElements[j]) {
                            isExplored[j] = true;
                        }
                    }
                }
                bottom *= factorial(sum);

            }
        }


        if(!dataKnown || !betaKnown){
            sum2=0;
        KPlus=0;
        for (int i = 0; i <data.getColumnDimension() ; i++) {
          if(containsNonZeroElements[i]) {
              KPlus++;
              for (int j = 0; j < data.getRowDimension(); j++) {
                  rowCount[i] += data.getParameterValue(j, i);
              }
              sum2+=Beta.logBeta(rowCount[i], data.getRowDimension() + beta.getParameterValue(0) - rowCount[i]);
          }
        }
        }
        double p1=KPlus*Math.log(alpha.getParameterValue(0)*beta.getParameterValue(0)/bottom);
        double p2=-alpha.getParameterValue(0)*H();
        double p3=sum2;
        betaKnown=true;
        dataKnown=true;
        return p1+p2+p3;
    }

    @Override
    public void makeDirty() {
        betaKnown=false;
        dataKnown=false;

    }

    boolean likelihoodKnown;
    boolean storedLikelihoodKnown;
    double logLikelihood;
    double storedLogLikelihood;
    boolean betaKnown=false;
    boolean dataKnown=false;
    boolean storedDataKnown;
    boolean storedBetaKnown;
    int[] rowCount;
    int[] storedRowCount;
    int KPlus;
    int storedKPlus;
    boolean[] containsNonZeroElements;
    boolean[] storedContainsNonZeroElements;
    double H;
    double storedH;
    int bottom;
    int storedBottom;
    double sum2;
    double storedSum2;

    MatrixParameter data;
    Parameter alpha;
    Parameter beta;
}
