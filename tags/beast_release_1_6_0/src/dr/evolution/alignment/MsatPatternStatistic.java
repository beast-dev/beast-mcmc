package dr.evolution.alignment;


/**
 * @author Chieh-Hsi
 */
public class MsatPatternStatistic {
    Patterns msatPattern;
    double msatLengthVar;
    String mode;


    public MsatPatternStatistic(Patterns msatPattern){
        this(msatPattern, "variance");
    }

    public MsatPatternStatistic(Patterns msatPattern, String mode){
        this.msatPattern = msatPattern;
        this.msatLengthVar = computeMsatLengthVariance();
        this.mode = mode;
    }

    public double computeMsatLengthVariance(){
        double var = 0.0;
        int[] msatPat = msatPattern.getPattern(0);
        double mean = 0.0;
        for(int i = 0; i < msatPat.length; i++){
            mean += msatPat[i];
        }
        mean = mean/msatPat.length;
        for(int i = 0; i < msatPat.length; i++){
            var+=(msatPat[i] - mean)*(msatPat[i] - mean);
        }
        var = var/msatPat.length;
        System.out.println(2*var);
        return var;

    }

    public String toString(){
        if(mode.equals("thetaV")){
            return ""+2*msatLengthVar;
        }else{
            return ""+msatLengthVar;
        }

    }


}
