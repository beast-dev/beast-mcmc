package dr.evomodel.epidemiology.casetocase;

/**
 * Created with IntelliJ IDEA.
 * User: mhall
 * Date: 22/08/2013
 * Time: 14:18
 * To change this template use File | Settings | File Templates.
 */
public class BadPartitionException extends RuntimeException {


    public BadPartitionException(String s){
        super(s);
    }

    public BadPartitionException(AbstractCase parentCase, AbstractCase childCase, double infectionTime){
        super("Suggesting that "+parentCase.getName()+" infected "+childCase.getName()+" at "+infectionTime+" which" +
                " is not permitted");
    }


}
