package dr.inference.model;


import org.apache.commons.math.special.Beta;

/**
 * @author Max Tolkoff
 */
public class IndianBuffetProcessPrior extends AbstractModelLikelihood {

    public IndianBuffetProcessPrior(Parameter alpha, Parameter beta, MatrixParameter data) {
        super(null);
        this.alpha=alpha;
        this.beta=beta;
        this.data=data;
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
        double sum=0;
        for (int i = 0; i <data.getRowDimension() ; i++) {
            sum+=beta.getParameterValue(0)/(beta.getParameterValue(0)+i);
        }
        return sum;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

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
        return calculateLogLikelihood();
    }

    private double calculateLogLikelihood(){
        int bottom=1;
        int sum;
        double sum2=0;
        boolean[] isExplored= new boolean[data.getColumnDimension()];
        boolean[] containsAll0=new boolean[data.getColumnDimension()];
        int[] rowCount=new int[data.getColumnDimension()];
        boolean same;
        for (int i = 0; i <data.getColumnDimension(); i++) {
            sum=1;
            for (int j = i+1; j <data.getColumnDimension() ; j++) {
                same=true;
                if(!isExplored[j]) {
                    for (int k = 0; k <data.getRowDimension(); k++){
                        if(data.getParameterValue(k,i)!=data.getParameterValue(k,j))
                            same=false;
                        if(data.getParameterValue(k,j)!=0){
                            containsAll0[j]=false;
                        }
//                        rowCount[j]+=data.getParameterValue(k,j);
                    }
                }
                if(same)
                {isExplored[j]=true;
                sum+=1;}

            }
            bottom*=factorial(sum);

        }
        int KPlus=0;
        for (int i = 0; i <data.getColumnDimension() ; i++) {
          if(!containsAll0[i]) {
              KPlus++;
              for (int j = 0; j < data.getRowDimension(); j++) {
                  rowCount[i] += data.getParameterValue(i, j);
              }
              sum2+=Beta.logBeta(rowCount[i], data.getRowDimension() + beta.getParameterValue(0) - rowCount[i]);
          }
        }
        double p1=KPlus*Math.log(alpha.getParameterValue(0)*beta.getParameterValue(0)/bottom);
        double p2=-alpha.getParameterValue(0)*H();
        double p3=sum2;
        return p1+p2+p3;
    }

    @Override
    public void makeDirty() {

    }

    MatrixParameter data;
    Parameter alpha;
    Parameter beta;
}
