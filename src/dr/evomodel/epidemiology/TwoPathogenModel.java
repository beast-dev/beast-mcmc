package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;

import java.util.List;

public class TwoPathogenModel extends CompartmentalModel {

    // did i put this in here?
    // private final DynamicalSystem rateModifiers;

    protected Parameter originTwo;
    protected Parameter originTimeNumSS;

    // To start, just coped code from SIRCompartmentalModel
    // Everything needs to be changed
    public TwoPathogenModel(
            List<Parameter> rateParams,
            List<Parameter> compartmentCounts,
            Parameter origin,
            Parameter originTwo,
            Parameter originTimeNumSS,
            int numReactionChannels,
            int numGridPoints,
            double cutOff) {
        super("Two-Pathogen CompartmentalModel");

        this.rateParameters = rateParams;
        for(int i = 0; i < rateParameters.size(); i++) {
            addVariable(rateParameters.get(i));
        }

        //this.rateModifiers = rateModifiers;
        //for(int i = 0; i < rateModifiers.size(); i++) {
        //    addVariable(rateModifiers.get(i));
        //}

        this.compartmentCounts = compartmentCounts;
        for(int i = 0; i < compartmentCounts.size(); i++) {
            addVariable(compartmentCounts.get(i));
        }
        this.origin = origin;
        addVariable(origin);
        this.originTwo = originTwo;
        addVariable(originTwo);
        this.originTimeNumSS = originTimeNumSS;
        addVariable(originTimeNumSS);
        this.numGridPoints = numGridPoints;
        this.cutOff = cutOff;
        this.numReactionChannels = numReactionChannels;
        this.numSpecies = compartmentCounts.size();
        this.vMatrix = getVMatrix();
        //this.originTimeNumS = originTimeNumS;
        //addVariable(originTimeNumS);
    }

    /*
    protected void setOriginTimeCompartmentCounts(int index){
        // initial SS value
        compartmentCounts.get(0).setParameterValue(index, 1000);
        // initial SI value
        compartmentCounts.get(1).setParameterValue(index, 1);
        // initial SC value
        compartmentCounts.get(2).setParameterValue(index, 0);
        // initial SR value
        compartmentCounts.get(3).setParameterValue(index, 0);
        // initial IS value
        compartmentCounts.get(4).setParameterValue(index, 0);
        // initial II value
        compartmentCounts.get(5).setParameterValue(index, 0);
        // initial IC value
        compartmentCounts.get(6).setParameterValue(index, 0);
        // initial IR value
        compartmentCounts.get(7).setParameterValue(index, 0);
        // initial CS value
        compartmentCounts.get(8).setParameterValue(index, 0);
        // initial CI value
        compartmentCounts.get(9).setParameterValue(index, 0);
        // initial CC value
        compartmentCounts.get(10).setParameterValue(index, 0);
        // initial CR value
        compartmentCounts.get(11).setParameterValue(index, 0);
        // initial RS value
        compartmentCounts.get(12).setParameterValue(index, 0);
        // initial RI value
        compartmentCounts.get(13).setParameterValue(index, 0);
        // initial RC value
        compartmentCounts.get(14).setParameterValue(index, 0);
        // initial RR value
        compartmentCounts.get(15).setParameterValue(index, 0);
    }

     */

    protected void setOriginTimeCompartmentCounts(int index){

        double origin1 = origin.getParameterValue(0);
        double origin2 = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);

        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 0);
        }

        // SS = originTimeSS
        compartmentCounts.get(0).setParameterValue(index, originTimeSS);

        if (origin1 > origin2) {
            // pathogen 1 is older, start in IS
            compartmentCounts.get(4).setParameterValue(index, 1);

        } else if (origin2 > origin1) {
            // pathogen 2 is older, start in SI
            compartmentCounts.get(1).setParameterValue(index, 1);

        } else {
            // origins equal
            // choose whichever convention you want
            compartmentCounts.get(1).setParameterValue(index, 1); // SI
            compartmentCounts.get(4).setParameterValue(index, 1); // IS
            compartmentCounts.get(0).setParameterValue(index, originTimeSS - 1); // starting with two sick
        }
    }

    /*
    protected void setDefaultCompartmentCounts(int index){
        // default S value
        compartmentCounts.get(0).setParameterValue(index, 999);
        // default I value
        compartmentCounts.get(1).setParameterValue(index, 0);
        // default R value
        compartmentCounts.get(2).setParameterValue(index, 0);
    }

     */

    protected void setDefaultCompartmentCounts(int index){
        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 0);
        }

        // default SS value
        compartmentCounts.get(0).setParameterValue(index, 999);


    }

    protected int[] getHighestOrdersOfReactions(){
        int[] gVec = new int[numSpecies];
        // g_i's determined by using the highest order of reaction of species i
        // g_SS = 2;
        gVec[0] = 2;
        // g_SI = 2;
        gVec[1] = 2;
        // g_SC = 2;
        gVec[2] = 2;
        // g_SR = 2;
        gVec[3] = 2;
        // g_IS = 2;
        gVec[4] = 2;
        // g_II = 2;
        gVec[5] = 2;
        // g_IC = 2;
        gVec[6] = 2;
        // g_IR = 2;
        gVec[7] = 2;
        // g_CS = 2;
        gVec[8] = 2;
        // g_CI = 2;
        gVec[9] = 2;
        // g_CC = 1;
        gVec[10] = 1;
        // g_CR = 1;
        gVec[11] = 1;
        // g_RS = 2;
        gVec[12] = 2;
        // g_RI = 2;
        gVec[13] = 2;
        // g_RC = 1;
        gVec[14] = 1;
        // g_RR = 1;
        gVec[15] = 1;
        return gVec;
    }

    // v matrix describes how count vector changes with reaction
    // row corresponds to species/particle type, column corresponds to reaction channel
    // columns are rxns, rows are compartment counts, there are 56 reactions
    // I am using the same reaction order as the reaction rates
    protected int[][] getVMatrix() {
        int[][] v = new int[][]{
                // SS
                {-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // SI
                { 0,  0,  0,  0,  1,  1,  1,  1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},

                // SC
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},

                // SR
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0},

                // IS
                { 1,  1,  1,  1,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // II
                { 0,  0,  0,  0,  0,  0,  0,  0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // IC
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // IR
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // CS
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // CI
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // CC
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // CR
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

                // RS
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1},

                // RI
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1, 0, 0, 0, 0},

                // RC
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 0, 0},

                // RR
                { 0,  0,  0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1}
        };

        return v;
    }

    protected double[] getReactionIntensities(double[] currentCounts){
        double[] rVec = new double[numReactionChannels];
        // currentCounts[0] has number of SS
        // currentCounts[1] has number of SI
        // currentCounts[2] has number of SC
        // currentCounts[3] has number of SR
        // currentCounts[4] has number of IS
        // currentCounts[5] has number of II
        // currentCounts[6] has number of IC
        // currentCounts[7] has number of IR
        // currentCounts[8] has number of CS
        // currentCounts[9] has number of CI
        // currentCounts[10] has number of CC
        // currentCounts[11] has number of CR
        // currentCounts[12] has number of RS
        // currentCounts[13] has number of RI
        // currentCounts[14] has number of RC
        // currentCounts[15] has number of RR
        //
        // transmissionRate[0] is transmission rate of pathogen 1
        // transmissionRate[1] is transmission rate of pathogen 2
        // same logic for recoveryRate, crossProtection and resusceptibilityRate
        // define crossProtection as I -> C and recoveryRate as C -> R and resusceptibilityRate as R _> S
        // look up how Chen et al defined each rate
        //
        // should the rates be rateParameters?
        // rateParameters: 0, 1 are transmissionRate 1 and 2; 2,3 are crossProtection 1 and 2; 4, 5 are recoveryRate 1 and 2
        //
        // rateModifiers[0] = alpha, modifies transmission during infection
        // rateModifiers[1] = chi, modifies transmission during cross protection
        // rateModifiers[2] = sigma, modifies recovery to cross protection during dual infection
        // one dimensional parameters
        // add crossProtection to parser and xml and 2 Two PathogenModel class
        // add rate modifiers as one dimensional parameters
        double transmissionRateOne = rateParameters.get(0).getParameterValue(0);
        double transmissionRateTwo = rateParameters.get(5).getParameterValue(0);
        double recoveryRateOne = rateParameters.get(1).getParameterValue(0);
        double recoveryRateTwo = rateParameters.get(6).getParameterValue(0);
        double crossProtectionRateOne = rateParameters.get(2).getParameterValue(0);
        double crossProtectionRateTwo = rateParameters.get(7).getParameterValue(0);
        double resusRateOne = rateParameters.get(4).getParameterValue(0);
        double resusRateTwo = rateParameters.get(9).getParameterValue(0);
        double infectionRateModulationI = rateParameters.get(10).getParameterValue(0);
        double infectionRateModulationC = rateParameters.get(11).getParameterValue(0);
        double recoveryRateModulation = rateParameters.get(12).getParameterValue(0);
        // SS + IS -> 2IS
        rVec[0] = transmissionRateOne*currentCounts[0]*currentCounts[4];
        // SS + II -> IS + II
        rVec[1] = transmissionRateOne*currentCounts[0]*currentCounts[5];
        // SS + IC -> IS + IC
        rVec[2] = transmissionRateOne*currentCounts[0]*currentCounts[6];
        // SS + IR -> IS + IR
        rVec[3] = transmissionRateOne*currentCounts[0]*currentCounts[7];
        // SS + SI -> 2SI
        rVec[4] = transmissionRateTwo*currentCounts[0]*currentCounts[1];
        // SS + II -> SI + II
        rVec[5] = transmissionRateTwo*currentCounts[0]*currentCounts[5];
        // SS + CI -> SI + CI
        rVec[6] = transmissionRateTwo*currentCounts[0]*currentCounts[9];
        // SS + RI -> SI + RI
        rVec[7] = transmissionRateTwo*currentCounts[0]*currentCounts[13];
        // SI + IS -> II + IS
        rVec[8] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[4];
        // SI + II -> 2II
        rVec[9] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[5];
        // SI + IC -> II + IC
        rVec[10] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[6];
        // SI + IR -> II + IR
        rVec[11] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[7];
        // SI -> SC
        rVec[12] = crossProtectionRateTwo*currentCounts[1];
        // SC + IS -> IC + IS
        rVec[13] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[4];
        // SC + II -> IC + II
        rVec[14] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[5];
        // SC + IC -> 2IC
        rVec[15] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[6];
        // SC + IR -> IC + IR
        rVec[16] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[7];
        // SC -> SR
        rVec[17] = recoveryRateTwo*currentCounts[2];
        // SR + IS -> IR + IS
        rVec[18] = transmissionRateOne*currentCounts[3]*currentCounts[4];
        // SR + II -> IR + II
        rVec[19] = transmissionRateOne*currentCounts[3]*currentCounts[5];
        // SR + IC -> IR + IC
        rVec[20] = transmissionRateOne*currentCounts[3]*currentCounts[6];
        // SR + IR -> 2IR
        rVec[21] = transmissionRateOne*currentCounts[3]*currentCounts[7];
        // SR -> SS
        rVec[22] = resusRateTwo*currentCounts[3];
        // IS -> CS
        rVec[23] = crossProtectionRateOne*currentCounts[4];
        // IS + SI -> II + SI
        rVec[24] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[1];
        // IS + II -> 2II
        rVec[25] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[5];
        // IS + CI -> II + CI
        rVec[26] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[9];
        // IS + RI -> II + RI
        rVec[27] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[13];
        // II -> CI
        rVec[28] = recoveryRateModulation*crossProtectionRateOne*currentCounts[5];
        // II -> IC
        rVec[29] = recoveryRateModulation*crossProtectionRateTwo*currentCounts[5];
        // IC -> CC
        rVec[30] = crossProtectionRateOne*currentCounts[5];
        // IC -> IR
        rVec[31] = recoveryRateTwo*currentCounts[5];
        // IR -> CR
        rVec[32] = crossProtectionRateOne*currentCounts[7];
        // IR -> IS
        rVec[33] = resusRateTwo*currentCounts[7];
        // CS -> RS
        rVec[34] = recoveryRateOne*currentCounts[8];
        // CS + SI -> CI + SI
        rVec[35] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[1];
        // CS + II -> CI + II
        rVec[36] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[5];
        // CS + CI -> 2CI
        rVec[37] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[9];
        // CS + RI -> CI + RI
        rVec[38] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[13];
        // CI -> RI
        rVec[39] = recoveryRateOne*currentCounts[9];
        // CI -> CC
        rVec[40] = crossProtectionRateTwo*currentCounts[9];
        // CC -> RC
        rVec[41] = recoveryRateOne*currentCounts[10];
        // CC -> CR
        rVec[42] = recoveryRateTwo*currentCounts[10];
        // CR -> RR
        rVec[43] = recoveryRateOne*currentCounts[11];
        // CR -> CS
        rVec[44] = resusRateTwo*currentCounts[11];
        // RS -> SS
        rVec[45] = resusRateOne*currentCounts[12];
        // RS + SI -> RI + SI
        rVec[46] = transmissionRateTwo*currentCounts[12]*currentCounts[1];
        // RS + II -> RI + II
        rVec[47] = transmissionRateTwo*currentCounts[12]*currentCounts[5];
        // RS + CI -> RI + CI
        rVec[48] = transmissionRateTwo*currentCounts[12]*currentCounts[9];
        // RS + RI -> 2RI
        rVec[49] = transmissionRateTwo*currentCounts[12]*currentCounts[13];
        // RI -> SI
        rVec[50] = resusRateOne*currentCounts[13];
        // RI -> RC
        rVec[51] = crossProtectionRateTwo*currentCounts[13];
        // RC -> SC
        rVec[52] = resusRateOne*currentCounts[14];
        // RC -> RR
        rVec[53] = recoveryRateTwo*currentCounts[14];
        // RR -> SR
        rVec[54] = resusRateOne*currentCounts[15];
        // RR -> RS
        rVec[55] = resusRateTwo*currentCounts[15];
        return rVec;
    }


    protected double[] getTimeDerivatives(double[] currentCounts){
        double[] returnVec = new double[numReactionChannels];
        double SS = currentCounts[0];
        double SI = currentCounts[1];
        double SC = currentCounts[2];
        double SR = currentCounts[3];
        double IS = currentCounts[4];
        double II = currentCounts[5];
        double IC = currentCounts[6];
        double IR = currentCounts[7];
        double CS = currentCounts[8];
        double CI = currentCounts[9];
        double CC = currentCounts[10];
        double CR = currentCounts[11];
        double RS = currentCounts[12];
        double RI = currentCounts[13];
        double RC = currentCounts[14];
        double RR = currentCounts[15];
        double transRateOne = rateParameters.get(0).getParameterValue(0);
        double recovRateOne = rateParameters.get(1).getParameterValue(0);
        double crossProtectionOne = rateParameters.get(2).getParameterValue(0);
        double resusRateOne = rateParameters.get(4).getParameterValue(0);
        double transRateTwo = rateParameters.get(5).getParameterValue(0);
        double recovRateTwo = rateParameters.get(6).getParameterValue(0);
        double crossProtectionTwo = rateParameters.get(7).getParameterValue(0);
        double resusRateTwo = rateParameters.get(9).getParameterValue(0);
        double alpha = rateParameters.get(10).getParameterValue(0);
        double chi = rateParameters.get(11).getParameterValue(0);
        double sigma = rateParameters.get(12).getParameterValue(0);
        // SS
        returnVec[0] = -transRateOne*SS*(IS + II + IR + IC) - transRateTwo*SS*(SI + II + CI + RI) + resusRateOne*RS + resusRateTwo*SR;
        // SI
        returnVec[1] = -alpha*transRateOne*SI*(IS + II + IC + IR) + resusRateOne*RI + transRateTwo*SS*(SI + II + CI + RI) - crossProtectionTwo*SI;
        // SC
        returnVec[2] = -chi*transRateOne*SC*(IS + II + IC + IR) + resusRateOne*RC + crossProtectionTwo*SI - recovRateTwo*SC;
        // SR
        returnVec[3] = -transRateOne*SR*(IS + II + IC + IR) + resusRateOne*RR - resusRateTwo*SR + recovRateTwo;
        // IS
        returnVec[4] = -alpha*transRateTwo*IS*(SI + II + CI + RI) + resusRateTwo*IR + transRateOne*SS*(IS + II + IC + IR) - crossProtectionOne;
        // II
        returnVec[5] = alpha*transRateOne*SI*(IS + II + CI + RI) - sigma*crossProtectionOne*II + alpha*transRateTwo*IS*(SI + II + CI + RI) - sigma*crossProtectionTwo*II;
        // IC
        returnVec[6] = chi*transRateTwo*SC*(IS + II + IC + IR) - crossProtectionOne*IC + sigma*crossProtectionTwo*II - recovRateTwo*IC;
        // IR
        returnVec[7] = transRateOne*SR*(IS + II + IC + IR) - crossProtectionOne*IR + recovRateTwo*IC - resusRateTwo*IR;
        // CS
        returnVec[8] = -chi*transRateTwo*CS*(SI + II + CI + RI) + resusRateTwo*CR + crossProtectionOne*IS - recovRateOne*CS;
        // CI
        returnVec[9] = chi*transRateTwo*CS*(SI + II + CI + RI) - crossProtectionTwo*CI + sigma*crossProtectionOne*II - recovRateOne*CI;
        // CC
        returnVec[10] = crossProtectionOne*IC - recovRateOne*CC + crossProtectionTwo*CI - recovRateTwo*CC;
        // CR
        returnVec[11] = crossProtectionOne*IR - recovRateOne*CR + recovRateTwo*CC - resusRateTwo*CR;
        // RS
        returnVec[12] = -transRateTwo*RS*(SI + II + CI + RI) + resusRateTwo*RR - resusRateOne*RS + recovRateOne*CS;
        // RI
        returnVec[13] = transRateTwo*RS*(SI + II + CI + RI) - crossProtectionTwo*RI + recovRateOne*CI - resusRateOne*RI;
        // RC
        returnVec[14] = crossProtectionTwo*RI - recovRateTwo*RC + recovRateOne*CC - resusRateOne*RC;
        // RR
        returnVec[15] = recovRateOne*CR - resusRateOne*RR + recovRateTwo*RC - resusRateTwo*RR;
        return returnVec;
    }


    protected double[] getSALPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau){
        double[] returnVal = new double[numReactionChannels];
        // for standard tau leaping
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = reactionInt[r]*tau;
        }
        // for SAL algorithm, need to add extra terms to account for linear change in intensity
        double[] timeDerivatives = getTimeDerivatives(currentCounts);
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = returnVal[r] + timeDerivatives[r]*tau*tau*0.5;
        }
        return returnVal;
    }
    // getSALPoissonIntensities should have dimension 56


    // countsNew has each rxn
    protected double[] getUpdatedCompartmentCounts(double[] currentCounts, double[] countsNew){
        double[] updatedCounts = new double[numSpecies];
        // SS(t+tau) = SS(t) - rxn0 - rxn1 - rxn2 - rxn3 - rxn4 - rxn5 - rxn6 - rxn7 + rxn22 + rxn45
        updatedCounts[0] = currentCounts[0] - countsNew[0] - countsNew[1] - countsNew[2] - countsNew[3] - countsNew[4] - countsNew[5] - countsNew[6] - countsNew[7] + countsNew[22] + countsNew[45];
        // SI(t+tau) = SI(t) - rxn8 - rxn9 - rxn10 - rxn11 - rxn12 + rxn4 + rxn50
        updatedCounts[1] = currentCounts[1] - countsNew[8] - countsNew[9] - countsNew[10] - countsNew[11] - countsNew[12] + countsNew[4] + countsNew[50];
        // SC(t+tau) = SC(t) - rxn13 - rxn14 - rxn15 - rxn16 - rxn17 + rxn12 + rxn52
        updatedCounts[2] = currentCounts[2] - countsNew[13] - countsNew[14] - countsNew[15] - countsNew[16] - countsNew[17] + countsNew[12] + countsNew[52];
        // SR(t+tau) = SR(t) - rxn18 - rxn19 - rxn20 -  rxn21 - rxn22 + rxn17 + rxn54
        updatedCounts[3] = currentCounts[3] - countsNew[18] - countsNew[19] - countsNew[20] - countsNew[21] - countsNew[22] + countsNew[17] + countsNew[54];
        // IS(t+tau) = IS(t) - rxn23 - rxn24 - rxn25 - rxn26 - rxn27 + rxn0 + rxn1 + rxn2 + rxn3 + rxn33
        updatedCounts[4] = currentCounts[4] - countsNew[23] - countsNew[24] - countsNew[25] - countsNew[26] - countsNew[27] + countsNew[0] + countsNew[1] + countsNew[2] + countsNew[3] + countsNew[33];
        // II(t+tau) = II(t) - rxn28 - rxn29 + rxn8 + rxn9 + rxn10 + rxn11 + rxn24 + rxn25 + rxn26 + rxn27
        updatedCounts[5] = currentCounts[5] - countsNew[28] - countsNew[29] + countsNew[8] + countsNew[9] + countsNew[10] + countsNew[11] + countsNew[24] + countsNew[25] + countsNew[26] + countsNew[27];
        // IC(t+tau) = IC(t) - rxn30 - rxn31 + rxn13 + rxn14 + rxn15 + rxn16 + rxn29
        updatedCounts[6] = currentCounts[6] - countsNew[30] - countsNew[31] + countsNew[13] + countsNew[14] + countsNew[15] + countsNew[16] + countsNew[29];
        // IR(t+tau) = IR(t) - rxn32 - rxn33 + rxn18 + rxn19 + rxn20 + rxn21 + rxn31
        updatedCounts[7] = currentCounts[7] - countsNew[32] - countsNew[33] + countsNew[18] + countsNew[19] + countsNew[20] + countsNew[21] + countsNew[31];
        // CS(t+tau) = CS(t) - rxn34 - rxn35 - rxn36 - rxn37 - rxn38 + rxn23 + rxn44
        updatedCounts[8] = currentCounts[8] - countsNew[34] - countsNew[35] - countsNew[36] - countsNew[37] - countsNew[38] + countsNew[23] + countsNew[44];
        // CI(t+tau) = CI(t) - rxn39 - rxn40 + rxn28 + rxn35 + rxn36 + rxn37 + rxn38
        updatedCounts[9] = currentCounts[9] - countsNew[39] - countsNew[40] + countsNew[28] + countsNew[35] + countsNew[36] + countsNew[37] + countsNew[38];
        // CC(t+tau) = CC(t) - rxn41 - rxn42 + rxn30 + rxn40
        updatedCounts[10] = currentCounts[10] - countsNew[41] - countsNew[42] + countsNew[30] + countsNew[40];
        // CR(t+tau) = CR(t) - rxn43 - rxn44 + rxn32 + rxn42
        updatedCounts[11] = currentCounts[11] - countsNew[43] - countsNew[44] + countsNew[32] + countsNew[42];
        // RS(t+tau) = RS(t) - rxn45 - rxn46 - rxn47 - rxn48 - rxn49 + rxn34 + rxn55
        updatedCounts[12] = currentCounts[12] - countsNew[45] - countsNew[46] - countsNew[47] - countsNew[48] - countsNew[49] + countsNew[34] + countsNew[55];
        // RI(t+tau) = RI(t) - rxn50 - rxn51 + rxn39 + rxn46 + rxn47 + rxn48 + rxn49
        updatedCounts[13] = currentCounts[13] - countsNew[50] - countsNew[51] + countsNew[39] + countsNew[46] + countsNew[47] + countsNew[48] + countsNew[49];
        // RC(t+tau) = RC(t) - rxn52 - rxn53 + rxn41 + rxn51
        updatedCounts[14] = currentCounts[14] - countsNew[52] - countsNew[53] + countsNew[41] + countsNew[51];
        // RR(t+tau) = RR(t) - rxn54 - rxn55 + rxn43 + rxn53
        updatedCounts[15] = currentCounts[15] - countsNew[54] - countsNew[55] + countsNew[43] + countsNew[53];
        return updatedCounts;
    }

    // add step to introduce 2nd pathogen
    private boolean secondPathogenIntroduced = false;

    // used AI for help debugging this. I added introduceSecondPathogen to compartmentalModel too but it doesn't do anything unless it is a 2 path model. I don't know if this is the right way to do it.

    @Override
    public double[] introduceSecondPathogen(
            double previousTime,
            double currentTime,
            double[] currentCounts) {

        double origin1 = origin.getParameterValue(0);
        double origin2 = originTwo.getParameterValue(0);

        double youngestOrigin = Math.min(origin1, origin2);

        if (!secondPathogenIntroduced &&
                // this should ensure that it only happens once
                previousTime > youngestOrigin &&
                currentTime <= youngestOrigin) {

            currentCounts[0] = currentCounts[0] - 1; // SS

            if (origin1 < origin2) {
                // pathogen 1 is younger, introduce pathogen 1
                currentCounts[4] = currentCounts[4] + 1; // IS
            } else {
                // pathogen 2 is younger, introduce pathogen 2
                currentCounts[1] = currentCounts[1] + 1; // SI
            }

            secondPathogenIntroduced = true;
            System.out.println("Both pathogens now active");
        }

        return currentCounts;
    }



}
