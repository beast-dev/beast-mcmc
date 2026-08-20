package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;

import java.util.List;

public class TwoPathogenModel extends CompartmentalModel {

    protected Parameter originTimeNumSS;
    protected double originTimeNumSI;
    protected double originTimeNumIS;
    protected double totalPopSize;
    protected boolean seasonalModel;
    protected double seasonalPeriod;
    protected double seasonalAmpOne;
    protected double seasonalAmpTwo;
    protected double seasonalPeakDayOne;
    protected double seasonalPeakDayTwo;
    // keep track of introduction of "younger" pathogen
    private boolean secondPathogenIntroduced;

    public TwoPathogenModel(
            List<Parameter> rateParams,
            List<Parameter> compartmentCounts,
            Parameter originOne,
            Parameter originTwo,
            Parameter originTimeNumSS,
            double originTimeNumSI,
            double originTimeNumIS,
            int numReactionChannels,
            int numGridPoints,
            double cutOff,
            boolean seasonalModel,
            Parameter seasonalPeriod,
            Parameter seasonalAmpOne,
            Parameter seasonalPeakDayOne,
            Parameter seasonalAmpTwo,
            Parameter seasonalPeakDayTwo) {
        super("Two-Pathogen CompartmentalModel");

        this.rateParameters = rateParams;
        for(int i = 0; i < rateParameters.size(); i++) {
            addVariable(rateParameters.get(i));
        }
        this.compartmentCounts = compartmentCounts;
        for(int i = 0; i < compartmentCounts.size(); i++) {
            addVariable(compartmentCounts.get(i));
        }
        this.originOne = originOne;
        addVariable(originOne);
        this.originTwo = originTwo;
        addVariable(originTwo);
        this.originTimeNumSS = originTimeNumSS;
        addVariable(originTimeNumSS);
        this.originTimeNumSI = originTimeNumSI;
        this.originTimeNumIS = originTimeNumIS;
        // changed later if necessary
        this.totalPopSize = originTimeNumSS.getParameterValue(0) + 1;
        this.numGridPoints = numGridPoints;
        this.cutOff = cutOff;
        this.numReactionChannels = numReactionChannels;
        this.numSpecies = compartmentCounts.size();
        this.vMatrix = getVMatrix();

        this.seasonalModel = seasonalModel;
        if(seasonalModel) {
            this.seasonalPeriod = seasonalPeriod.getParameterValue(0);
            this.seasonalAmpOne = seasonalAmpOne.getParameterValue(0);
            this.seasonalAmpTwo = seasonalAmpTwo.getParameterValue(0);
            this.seasonalAmpTwo = seasonalAmpTwo.getParameterValue(0);
            this.seasonalPeakDayOne = seasonalAmpTwo.getParameterValue(0);
            this.seasonalPeakDayTwo = seasonalAmpTwo.getParameterValue(0);
        }
    }


    protected void setOriginTimeCompartmentCounts(int index){
        secondPathogenIntroduced = false;

        double origOne = originOne.getParameterValue(0);
        double origTwo = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);
        double originTimeSI = originTimeNumSI;
        double originTimeIS = originTimeNumIS;

        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 0);
        }

        // SS = originTimeSS
        compartmentCounts.get(0).setParameterValue(index, originTimeSS);

        if (origOne > origTwo) {
            // pathogen 1 is older, start in IS
            compartmentCounts.get(4).setParameterValue(index, originTimeIS);
            // total compartment counts at origin time should be originTimeSS + originTimeIS
            totalPopSize = originTimeSS + originTimeIS;
        } else if (origTwo > origOne) {
            // pathogen 2 is older, start in SI
            compartmentCounts.get(1).setParameterValue(index, originTimeSI);
            // total compartment counts at origin time should be originTimeSS + originTimeSI
            totalPopSize = originTimeSS + originTimeSI;
        } else {
            // no need to "introduce" second pathogen while doing forward time simulation
            secondPathogenIntroduced = true;
            // origins equal
            compartmentCounts.get(1).setParameterValue(index, originTimeSI); // SI
            compartmentCounts.get(4).setParameterValue(index, originTimeIS); // IS
            // total compartment counts at origin time should be originTimeSS + originTimeIS + originTimeSI
            totalPopSize = originTimeSS + originTimeIS + originTimeSI;
        }
    }

    /*
    // DEBUG AND SIMULATION ONLY - SET ALL COUNTS TO 50 INITIALLY
    protected void setOriginTimeCompartmentCounts(int index){

        double origOne = originOne.getParameterValue(0);
        double origTwo = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);

        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 50);
        }

        // total compartment counts should be originTimeSS + 1 for the infected individual?
        // SS = originTimeSS
        compartmentCounts.get(0).setParameterValue(index, originTimeSS);

        if (origOne > origTwo) {
            // pathogen 1 is older, start in IS
            compartmentCounts.get(4).setParameterValue(index, 50);
        } else if (origTwo > origOne) {
            // pathogen 2 is older, start in SI
            compartmentCounts.get(1).setParameterValue(index, 50);
        } else {
            // no need to "introduce" second pathogen while doing forward time simulation
            secondPathogenIntroduced = true;
            // origins equal
            // choose whichever convention you want
            compartmentCounts.get(1).setParameterValue(index, 50); // SI
            compartmentCounts.get(4).setParameterValue(index, 50); // IS
            compartmentCounts.get(0).setParameterValue(index, originTimeSS - 1); // starting with two sick
        }
    }
    */

    protected void setDefaultCompartmentCounts(int index){
        double origOne = originOne.getParameterValue(0);
        double origTwo = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);
        double originTimeSI = originTimeNumSI;
        double originTimeIS = originTimeNumIS;

        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 0);
        }

        if (origOne > origTwo) {
            // pathogen 1 is older, start in IS
            // total compartment counts at origin time should be originTimeSS + originTimeIS
            // default SS value should be same
            compartmentCounts.get(0).setParameterValue(index, originTimeSS + originTimeIS);
        } else if (origTwo > origOne) {
            // pathogen 2 is older, start in SI
            // total compartment counts at origin time should be originTimeSS + originTimeSI
            // default SS value should be same
            compartmentCounts.get(0).setParameterValue(index, originTimeSS + originTimeSI);
        } else {
            // no need to "introduce" second pathogen while doing forward time simulation
            secondPathogenIntroduced = true;
            // origins equal
            // total compartment counts at origin time should be originTimeSS + originTimeIS + originTimeSI
            // default value should be same
            compartmentCounts.get(0).setParameterValue(index, originTimeSS + originTimeIS + originTimeSI);
        }
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
        int[][] v = new int[numSpecies][numReactionChannels];

        // SS
        v[0][0] = -1; v[0][1] = -1; v[0][2] = -1; v[0][3] = -1; v[0][16] = -1; v[0][17] = -1; v[0][18] = -1;
        v[0][19] = -1; v[0][48] = 1; v[0][52] = 1; v[0][56] = -1; v[0][60] = -1; v[0][64] = 1; v[0][65] = 1;
        v[0][66] = 1; v[0][67] = 1; v[0][68] = 1; v[0][69] = 1; v[0][70] = 1; v[0][71] = 1; v[0][72] = 1;
        v[0][73] = 1; v[0][74] = 1; v[0][75] = 1; v[0][76] = 1; v[0][77] = 1; v[0][78] = 1;

        // SI
        v[1][4] = -1; v[1][5] = -1; v[1][6] = -1; v[1][7] = -1; v[1][16] = 1; v[1][17] = 1; v[1][18] = 1;
        v[1][19] = 1; v[1][36] = -1; v[1][49] = 1; v[1][57] = -1; v[1][60] = 1; v[1][67] = -1;

        // SC
        v[2][8] = -1; v[2][9] = -1; v[2][10] = -1; v[2][11] = -1; v[2][36] = 1;
        v[2][44] = -1; v[2][50] = 1; v[2][58] = -1; v[2][71] = -1;

        // SR
        v[3][12] = -1; v[3][13] = -1; v[3][14] = -1; v[3][15] = -1; v[3][44] = 1;
        v[3][51] = 1; v[3][52] = -1; v[3][59] = -1; v[3][75] = -1;

        // IS
        v[4][0] = 1; v[4][1] = 1; v[4][2] = 1; v[4][3] = 1; v[4][20] = -1; v[4][21] = -1; v[4][22] = -1;
        v[4][23] = -1; v[4][32] = -1; v[4][53] = 1; v[4][56] = 1; v[4][61] = -1; v[4][64] = -1;

        // II
        v[5][4] = 1; v[5][5] = 1; v[5][6] = 1; v[5][7] = 1; v[5][20] = 1; v[5][21] = 1; v[5][22] = 1;
        v[5][23] = 1; v[5][33] = -1; v[5][37] = -1; v[5][57] = 1; v[5][61] = 1; v[5][68] = -1;

        // IC
        v[6][8] = 1; v[6][9] = 1; v[6][10] = 1; v[6][11] = 1; v[6][34] = -1; v[6][37] = 1; v[6][45] = -1;
        v[6][58] = 1; v[6][72] = -1;

        // IR
        v[7][12] = 1; v[7][13] = 1; v[7][14] = 1; v[7][15] = 1; v[7][35] = -1; v[7][45] = 1;
        v[7][53] = -1; v[7][59] = 1; v[7][76] = -1;

        // CS
        v[8][24] = -1; v[8][25] = -1; v[8][26] = -1; v[8][27] = -1; v[8][32] = 1; v[8][40] = -1;
        v[8][54] = 1; v[8][62] = -1; v[8][65] = -1;

        // CI
        v[9][24] = 1; v[9][25] = 1; v[9][26] = 1; v[9][27] = 1; v[9][33] = 1; v[9][38] = -1;
        v[9][41] = -1; v[9][62] = 1; v[9][69] = -1;

        // CC
        v[10][34] = 1; v[10][38] = 1; v[10][42] = -1; v[10][46] = -1; v[10][73] = -1;

        // CR
        v[11][35] = 1; v[11][43] = -1; v[11][46] = 1; v[11][54] = -1; v[11][77] = -1;

        // RS
        v[12][28] = -1; v[12][29] = -1; v[12][30] = -1; v[12][31] = -1; v[12][40] = 1; v[12][48] = -1;
        v[12][55] = 1; v[12][63] = -1; v[12][66] = -1;

        // RI
        v[13][28] = 1; v[13][29] = 1; v[13][30] = 1; v[13][31] = 1; v[13][39] = -1; v[13][41] = 1;
        v[13][49] = -1; v[13][63] = 1; v[13][70] = -1;

        // RC
        v[14][39] = 1; v[14][42] = 1; v[14][47] = -1; v[14][50] = -1; v[14][74] = -1;

        // RR
        v[15][43] = 1; v[15][47] = 1; v[15][51] = -1; v[15][55] = -1; v[15][78] = -1;

        /*
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
        */
        return v;
    }

    protected double[] getReactionIntensities(double[] currentCounts, double simTime){
        double[] rVec = new double[numReactionChannels];
        double numSS = currentCounts[0];
        double numSI = currentCounts[1];
        double numSC = currentCounts[2];
        double numSR = currentCounts[3];
        double numIS = currentCounts[4];
        double numII = currentCounts[5];
        double numIC = currentCounts[6];
        double numIR = currentCounts[7];
        double numCS = currentCounts[8];
        double numCI = currentCounts[9];
        double numCC = currentCounts[10];
        double numCR = currentCounts[11];
        double numRS = currentCounts[12];
        double numRI = currentCounts[13];
        double numRC = currentCounts[14];
        double numRR = currentCounts[15];

        //double transmissionRateOne = rateParameters.get(0).getParameterValue(0);
        double transmissionRateOne = getTransmissionRateOne(simTime);
        double moveToCRateOne = rateParameters.get(1).getParameterValue(0);
        double moveToRRateOne = rateParameters.get(2).getParameterValue(0);
        double resusRateOne = rateParameters.get(4).getParameterValue(0);
        //double transmissionRateTwo = rateParameters.get(5).getParameterValue(0);
        double transmissionRateTwo = getTransmissionRateTwo(simTime);
        double moveToCRateTwo = rateParameters.get(6).getParameterValue(0);
        double moveToRRateTwo = rateParameters.get(7).getParameterValue(0);
        double resusRateTwo = rateParameters.get(9).getParameterValue(0);
        double infectionRateModulationI = rateParameters.get(10).getParameterValue(0);
        double infectionRateModulationC = rateParameters.get(11).getParameterValue(0);
        double infectionRateModulationR = rateParameters.get(12).getParameterValue(0);
        double birthMortalityRate = rateParameters.get(13).getParameterValue(0);
        double importationRateOne = rateParameters.get(14).getParameterValue(0);
        double importationRateTwo = rateParameters.get(15).getParameterValue(0);

        // Infection Reactions

        // SS + IS -> 2IS
        rVec[0] = transmissionRateOne*numSS*numIS;
        // SS + II -> IS + II
        rVec[1] = transmissionRateOne*numSS*numII;
        // SS + IC -> IS + IC
        rVec[2] = transmissionRateOne*numSS*numIC;
        // SS + IR -> IS + IR
        rVec[3] = transmissionRateOne*numSS*numIR;

        // SI + IS -> II + IS
        rVec[4] = infectionRateModulationI*transmissionRateOne*numSI*numIS;
        // SI + II -> 2II
        rVec[5] = infectionRateModulationI*transmissionRateOne*numSI*numII;
        // SI + IC -> II + IC
        rVec[6] = infectionRateModulationI*transmissionRateOne*numSI*numIC;
        // SI + IR -> II + IR
        rVec[7] = infectionRateModulationI*transmissionRateOne*numSI*numIR;

        // SC + IS -> IC + IS
        rVec[8] = infectionRateModulationC*transmissionRateOne*numSC*numIS;
        // SC + II -> IC + II
        rVec[9] = infectionRateModulationC*transmissionRateOne*numSC*numII;
        // SC + IC -> 2IC
        rVec[10] = infectionRateModulationC*transmissionRateOne*numSC*numIC;
        // SC + IR -> IC + IR
        rVec[11] = infectionRateModulationC*transmissionRateOne*numSC*numIR;

        // SR + IS -> IR + IS
        rVec[12] = infectionRateModulationR*transmissionRateOne*numSR*numIS;
        // SR + II -> IR + II
        rVec[13] = infectionRateModulationR*transmissionRateOne*numSR*numII;
        // SR + IC -> IR + IC
        rVec[14] = infectionRateModulationR*transmissionRateOne*numSR*numIC;
        // SR + IR -> 2IR
        rVec[15] = infectionRateModulationR*transmissionRateOne*numSR*numIR;

        // SS + SI -> 2SI
        rVec[16] = transmissionRateTwo*numSS*numSI;
        // SS + II -> SI + II
        rVec[17] = transmissionRateTwo*numSS*numII;
        // SS + CI -> SI + CI
        rVec[18] = transmissionRateTwo*numSS*numCI;
        // SS + RI -> SI + RI
        rVec[19] = transmissionRateTwo*numSS*numRI;

        // IS + SI -> II + SI
        rVec[20] = infectionRateModulationI*transmissionRateTwo*numIS*numSI;
        // IS + II -> 2II
        rVec[21] = infectionRateModulationI*transmissionRateTwo*numIS*numII;
        // IS + CI -> II + CI
        rVec[22] = infectionRateModulationI*transmissionRateTwo*numIS*numCI;
        // IS + RI -> II + RI
        rVec[23] = infectionRateModulationI*transmissionRateTwo*numIS*numRI;

        // CS + SI -> CI + SI
        rVec[24] = infectionRateModulationC*transmissionRateTwo*numCS*numSI;
        // CS + II -> CI + II
        rVec[25] = infectionRateModulationC*transmissionRateTwo*numCS*numII;
        // CS + CI -> 2CI
        rVec[26] = infectionRateModulationC*transmissionRateTwo*numCS*numCI;
        // CS + RI -> CI + RI
        rVec[27] = infectionRateModulationC*transmissionRateTwo*numCS*numRI;

        // RS + SI -> RI + SI
        rVec[28] = infectionRateModulationR*transmissionRateTwo*numRS*numSI;
        // RS + II -> RI + II
        rVec[29] = infectionRateModulationR*transmissionRateTwo*numRS*numII;
        // RS + CI -> RI + CI
        rVec[30] = infectionRateModulationR*transmissionRateTwo*numRS*numCI;
        // RS + RI -> 2RI
        rVec[31] = infectionRateModulationR*transmissionRateTwo*numRS*numRI;

        // Recovery reactions

        // IS -> CS
        rVec[32] = moveToCRateOne*numIS;
        // II -> CI
        rVec[33] = moveToCRateOne*numII;
        // IC -> CC
        rVec[34] = moveToCRateOne*numIC;
        // IR -> CR
        rVec[35] = moveToCRateOne*numIR;
        // SI -> SC
        rVec[36] = moveToCRateTwo*numSI;
        // II -> IC
        rVec[37] = moveToCRateTwo*numII;
        // CI -> CC
        rVec[38] = moveToCRateTwo*numCI;
        // RI -> RC
        rVec[39] = moveToCRateTwo*numRI;

        // Loss of cross-immunity reactions

        // CS -> RS
        rVec[40] = moveToRRateOne*numCS;
        // CI -> RI
        rVec[41] = moveToRRateOne*numCI;
        // CC -> RC
        rVec[42] = moveToRRateOne*numCC;
        // CR -> RR
        rVec[43] = moveToRRateOne*numCR;
        // SC -> SR
        rVec[44] = moveToRRateTwo*numSC;
        // IC -> IR
        rVec[45] = moveToRRateTwo*numIC;
        // CC -> CR
        rVec[46] = moveToRRateTwo*numCC;
        // RC -> RR
        rVec[47] = moveToRRateTwo*numRC;

        // Resusceptibility reactions

        // RS -> SS
        rVec[48] = resusRateOne*numRS;
        // RI -> SI
        rVec[49] = resusRateOne*numRI;
        // RC -> SC
        rVec[50] = resusRateOne*numRC;
        // RR -> SR
        rVec[51] = resusRateOne*numRR;
        // SR -> SS
        rVec[52] = resusRateTwo*numSR;
        // IR -> IS
        rVec[53] = resusRateTwo*numIR;
        // CR -> CS
        rVec[54] = resusRateTwo*numCR;
        // RR -> RS
        rVec[55] = resusRateTwo*numRR;

        // Infection due to importation

        // SS -> IS
        rVec[56] = transmissionRateOne*importationRateOne*numSS;
        // SI -> II
        rVec[57] = infectionRateModulationI*transmissionRateOne*importationRateOne*numSI;
        // SC -> IC
        rVec[58] = infectionRateModulationC*transmissionRateOne*importationRateOne*numSC;
        // SR -> IR
        rVec[59] = infectionRateModulationR*transmissionRateOne*importationRateOne*numSR;
        // SS -> SI
        rVec[60] = transmissionRateTwo*importationRateTwo*numSS;
        // IS -> II
        rVec[61] = infectionRateModulationI*transmissionRateTwo*importationRateTwo*numIS;
        // CS -> CI
        rVec[62] = infectionRateModulationC*transmissionRateTwo*importationRateTwo*numCS;
        // RS -> RI
        rVec[63] = infectionRateModulationR*transmissionRateTwo*importationRateTwo*numRS;


        // Demographic replacement (birth-death reactions)
        // NOTE: we are omitting the trivial SS -> 0, 0 -> SS reaction, so there are 79 reactions in total
        // rather than the 80 that there should be in theory

        // IS -> 0, 0 -> SS
        rVec[64] = birthMortalityRate*numIS;
        // CS -> 0, 0 -> SS
        rVec[65] = birthMortalityRate*numCS;
        // RS -> 0, 0 -> SS
        rVec[66] = birthMortalityRate*numRS;
        // SI -> 0, 0 -> SS
        rVec[67] = birthMortalityRate*numSI;
        // II -> 0, 0 -> SS
        rVec[68] = birthMortalityRate*numII;
        // CI -> 0, 0 -> SS
        rVec[69] = birthMortalityRate*numCI;
        // RI -> 0, 0 -> SS
        rVec[70] = birthMortalityRate*numRI;
        // SC -> 0, 0 -> SS
        rVec[71] = birthMortalityRate*numSC;
        // IC -> 0, 0 -> SS
        rVec[72] = birthMortalityRate*numIC;
        // CC -> 0, 0 -> SS
        rVec[73] = birthMortalityRate*numCC;
        // RC -> 0, 0 -> SS
        rVec[74] = birthMortalityRate*numRC;
        // SR -> 0, 0 -> SS
        rVec[75] = birthMortalityRate*numSR;
        // IR -> 0, 0 -> SS
        rVec[76] = birthMortalityRate*numIR;
        // CR -> 0, 0 -> SS
        rVec[77] = birthMortalityRate*numCR;
        // RR -> 0, 0 -> SS
        rVec[78] = birthMortalityRate*numRR;

        /*
        // SS + IS -> 2IS
        rVec[0] = transmissionRateOne*numSS*numIS;
        // SS + II -> IS + II
        rVec[1] = transmissionRateOne*numSS*numII;
        // SS + IC -> IS + IC
        rVec[2] = transmissionRateOne*numSS*numIC;
        // SS + IR -> IS + IR
        rVec[3] = transmissionRateOne*numSS*numIR;
        // SS + SI -> 2SI
        rVec[4] = transmissionRateTwo*numSS*numSI;
        // SS + II -> SI + II
        rVec[5] = transmissionRateTwo*numSS*numII;
        // SS + CI -> SI + CI
        rVec[6] = transmissionRateTwo*numSS*numCI;
        // SS + RI -> SI + RI
        rVec[7] = transmissionRateTwo*numSS*numRI;
        // SI + IS -> II + IS
        rVec[8] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[4];
        // SI + II -> 2II
        rVec[9] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[5];
        // SI + IC -> II + IC
        rVec[10] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[6];
        // SI + IR -> II + IR
        rVec[11] = infectionRateModulationI*transmissionRateOne*currentCounts[1]*currentCounts[7];
        // SI -> SC
        rVec[12] = moveToCRateTwo *currentCounts[1];
        // SC + IS -> IC + IS
        rVec[13] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[4];
        // SC + II -> IC + II
        rVec[14] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[5];
        // SC + IC -> 2IC
        rVec[15] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[6];
        // SC + IR -> IC + IR
        rVec[16] = infectionRateModulationC*transmissionRateOne*currentCounts[2]*currentCounts[7];
        // SC -> SR
        rVec[17] = moveToRRateTwo *currentCounts[2];
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
        rVec[23] = moveToCRateOne*currentCounts[4];
        // IS + SI -> II + SI
        rVec[24] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[1];
        // IS + II -> 2II
        rVec[25] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[5];
        // IS + CI -> II + CI
        rVec[26] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[9];
        // IS + RI -> II + RI
        rVec[27] = infectionRateModulationI*transmissionRateTwo*currentCounts[4]*currentCounts[13];
        // II -> CI
        rVec[28] = infectionRateModulationR*moveToCRateOne*currentCounts[5];
        // II -> IC
        rVec[29] = infectionRateModulationR*moveToCRateTwo*currentCounts[5];
        // IC -> CC
        rVec[30] = moveToCRateOne*currentCounts[6];
        // IC -> IR
        rVec[31] = moveToRRateTwo *currentCounts[6];
        // IR -> CR
        rVec[32] = moveToCRateOne*currentCounts[7];
        // IR -> IS
        rVec[33] = resusRateTwo*currentCounts[7];
        // CS -> RS
        rVec[34] = moveToRRateOne *currentCounts[8];
        // CS + SI -> CI + SI
        rVec[35] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[1];
        // CS + II -> CI + II
        rVec[36] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[5];
        // CS + CI -> 2CI
        rVec[37] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[9];
        // CS + RI -> CI + RI
        rVec[38] = infectionRateModulationC*transmissionRateTwo*currentCounts[8]*currentCounts[13];
        // CI -> RI
        rVec[39] = moveToRRateOne*currentCounts[9];
        // CI -> CC
        rVec[40] = moveToCRateTwo*currentCounts[9];
        // CC -> RC
        rVec[41] = moveToRRateOne*currentCounts[10];
        // CC -> CR
        rVec[42] = moveToRRateTwo*currentCounts[10];
        // CR -> RR
        rVec[43] = moveToRRateOne*currentCounts[11];
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
        rVec[51] = moveToCRateTwo *currentCounts[13];
        // RC -> SC
        rVec[52] = resusRateOne*currentCounts[14];
        // RC -> RR
        rVec[53] = moveToRRateTwo *currentCounts[14];
        // RR -> SR
        rVec[54] = resusRateOne*currentCounts[15];
        // RR -> RS
        rVec[55] = resusRateTwo*currentCounts[15];
        */
        return rVec;
    }


    protected double[] getCompartmentDerivatives(double[] currentCounts, double simTime){

        double[] returnVec = new double[numSpecies];

        double numSS = currentCounts[0];
        double numSI = currentCounts[1];
        double numSC = currentCounts[2];
        double numSR = currentCounts[3];
        double numIS = currentCounts[4];
        double numII = currentCounts[5];
        double numIC = currentCounts[6];
        double numIR = currentCounts[7];
        double numCS = currentCounts[8];
        double numCI = currentCounts[9];
        double numCC = currentCounts[10];
        double numCR = currentCounts[11];
        double numRS = currentCounts[12];
        double numRI = currentCounts[13];
        double numRC = currentCounts[14];
        double numRR = currentCounts[15];

        //double transRateOne = rateParameters.get(0).getParameterValue(0);
        double moveToCRateOne = rateParameters.get(1).getParameterValue(0);
        double moveToRRateOne = rateParameters.get(2).getParameterValue(0);
        double resusRateOne = rateParameters.get(4).getParameterValue(0);
        //double transRateTwo = rateParameters.get(5).getParameterValue(0);
        double moveToCRateTwo = rateParameters.get(6).getParameterValue(0);
        double moveToRRateTwo = rateParameters.get(7).getParameterValue(0);
        double resusRateTwo = rateParameters.get(9).getParameterValue(0);

        double infectionRateModulationI = rateParameters.get(10).getParameterValue(0);
        double infectionRateModulationC = rateParameters.get(11).getParameterValue(0);
        double infectionRateModulationR = rateParameters.get(12).getParameterValue(0);

        double birthMortalityRate = rateParameters.get(13).getParameterValue(0);
        //double importationRateOne = rateParameters.get(14).getParameterValue(0);
        //double importationRateTwo = rateParameters.get(15).getParameterValue(0);

        double forceOfInfectionRateOne = getForceOfInfectionOne(currentCounts, simTime);
        double forceOfInfectionRateTwo = getForceOfInfectionTwo(currentCounts, simTime);

        // SS
        //returnVec[0] = -transRateOne*numSS*(numIS + numII + numIR + numIC)
        //        - transRateTwo*numSS*(numSI + numII + numCI + numRI)
        //        + resusRateOne*numRS + resusRateTwo*numSR;
        returnVec[0] = birthMortalityRate*totalPopSize - forceOfInfectionRateOne*numSS
                - forceOfInfectionRateTwo*numSS + resusRateOne*numRS + resusRateTwo*numSR
                - birthMortalityRate*numSS;
        // SI
        //returnVec[1] = -infectionRateModulationI*transRateOne*numSI*(numIS + numII + numIC + numIR)
        //        + resusRateOne*numRI
        //        + transRateTwo*numSS*(numSI + numII + numCI + numRI) - moveToCRateTwo*numSI;
        returnVec[1] = -infectionRateModulationI*forceOfInfectionRateOne*numSI + resusRateOne*numRI
                + forceOfInfectionRateTwo*numSS - moveToCRateTwo*numSI - birthMortalityRate*numSI;
        // SC
        //returnVec[2] = -infectionRateModulationC*transRateOne*numSC*(numIS + numII + numIC + numIR)
        //        + resusRateOne*numRC
        //        + moveToCRateTwo*numSI - moveToRRateTwo*numSC;
        returnVec[2] = -infectionRateModulationC*forceOfInfectionRateOne*numSC + resusRateOne*numRC
                + moveToCRateTwo*numSI - moveToRRateTwo*numSC - birthMortalityRate*numSC;
        // SR
        //returnVec[3] = -transRateOne*numSR*(numIS + numII + numIC + numIR)
        //        + resusRateOne*numRR - resusRateTwo*numSR
        //        + moveToRRateTwo;
        returnVec[3] = -infectionRateModulationR*forceOfInfectionRateOne*numSR
                + resusRateOne*numRR - resusRateTwo*numSR
                + moveToRRateTwo*numSC - birthMortalityRate*numSR;
        // IS
        //returnVec[4] = -infectionRateModulationI*transRateTwo*numIS*(numSI + numII + numCI + numRI)
        //        + resusRateTwo*numIR
        //        + transRateOne*numSS*(numIS + numII + numIC + numIR) - moveToCRateOne;
        returnVec[4] = -infectionRateModulationI*forceOfInfectionRateTwo*numIS
                + resusRateTwo*numIR + forceOfInfectionRateOne*numSS - moveToCRateOne*numIS
                - birthMortalityRate*numIS;
        // II
        //returnVec[5] = infectionRateModulationI*transRateOne*numSI*(numIS + numII + numCI + numRI)
        //        - infectionRateModulationR* moveToCRateOne*numII
        //        + infectionRateModulationI*transRateTwo*numIS*(numSI + numII + numCI + numRI)
        //        - infectionRateModulationR*moveToCRateTwo*numII;
        returnVec[5] = infectionRateModulationI*forceOfInfectionRateOne*numSI
                - moveToCRateOne*numII + infectionRateModulationI*forceOfInfectionRateTwo*numIS
                - moveToCRateTwo*numII - birthMortalityRate*numII;
        // IC
        //returnVec[6] = infectionRateModulationC*transRateTwo*numSC*(numIS + numII + numIC + numIR)
        //        - moveToCRateOne*numIC + infectionRateModulationR*moveToCRateTwo*numII - moveToRRateTwo*numIC;
        returnVec[6] = infectionRateModulationC*forceOfInfectionRateOne*numSC - moveToCRateOne*numIC
                + moveToCRateTwo*numII - moveToRRateTwo*numIC - birthMortalityRate*numIC;
        // IR
        //returnVec[7] = transRateOne*numSR*(numIS + numII + numIC + numIR)
        //        - moveToCRateOne*numIR + moveToRRateTwo*numIC - resusRateTwo*numIR;
        returnVec[7] = infectionRateModulationR*forceOfInfectionRateOne*numSR - moveToCRateOne*numIR
                + moveToRRateTwo*numIC - resusRateTwo*numIR - birthMortalityRate*numIR;
        // CS
        //returnVec[8] = -infectionRateModulationC*transRateTwo*numCS*(numSI + numII + numCI + numRI)
        //        + resusRateTwo*numCR + moveToCRateOne*numIS - moveToRRateOne*numCS;
        returnVec[8] = -infectionRateModulationC*forceOfInfectionRateTwo*numCS + resusRateTwo*numCR
                + moveToCRateOne*numIS - moveToRRateOne*numCS - birthMortalityRate*numCS;
        // CI
        //returnVec[9] = infectionRateModulationC*transRateTwo*numCS*(numSI + numII + numCI + numRI)
        //        - moveToCRateTwo*numCI + infectionRateModulationR*moveToCRateOne*numII - moveToRRateOne*numCI;
        returnVec[9] = infectionRateModulationC*forceOfInfectionRateTwo*numCS - moveToCRateTwo*numCI + moveToCRateOne*numII
                - moveToRRateOne*numCI - birthMortalityRate*numCI;
        // CC
        //returnVec[10] = moveToCRateOne*numIC - moveToRRateOne*numCC + moveToCRateTwo*numCI - moveToRRateTwo*numCC;
        returnVec[10] = moveToCRateOne*numIC - moveToRRateOne*numCC + moveToCRateTwo*numCI
                - moveToRRateTwo*numCC - birthMortalityRate*numCC;
        // CR
        //returnVec[11] = moveToCRateOne*numIR - moveToRRateOne*numCR + moveToRRateTwo*numCC - resusRateTwo*numCR;
        returnVec[11] = moveToCRateOne*numIR - moveToRRateOne*numCR + moveToRRateTwo*numCC
                - resusRateTwo*numCR - birthMortalityRate*numCR;
        // RS
        //returnVec[12] = -transRateTwo*numRS*(numSI + numII + numCI + numRI)
        //        + resusRateTwo*numRR - resusRateOne*numRS + moveToRRateOne*numCS;
        returnVec[12] = -infectionRateModulationR*forceOfInfectionRateTwo*numRS + resusRateTwo*numRR
                - resusRateOne*numRS + moveToRRateOne*numCS - birthMortalityRate*numRS;
        // RI
        //returnVec[13] = transRateTwo*numRS*(numSI + numII + numCI + numRI)
        //        - moveToCRateTwo*numRI + moveToRRateOne*numCI - resusRateOne*numRI;
        returnVec[13] = infectionRateModulationR*forceOfInfectionRateTwo*numRS - moveToCRateTwo*numRI
                + moveToRRateOne*numCI - resusRateOne*numRI - birthMortalityRate*numRI;
        // RC
        //returnVec[14] = moveToCRateTwo*numRI - moveToRRateTwo*numRC + moveToRRateOne*numCC - resusRateOne*numRC;
        returnVec[14] = moveToCRateTwo*numRI - moveToRRateTwo*numRC + moveToRRateOne*numCC
                - resusRateOne*numRC - birthMortalityRate*numRC;
        // RR
        //returnVec[15] = moveToRRateOne*numCR - resusRateOne*numRR + moveToRRateTwo*numRC - resusRateTwo*numRR;
        returnVec[15] = moveToRRateOne*numCR - resusRateOne*numRR + moveToRRateTwo*numRC
                - resusRateTwo*numRR - birthMortalityRate*numRR;
        return returnVec;
    }

    protected double[] getTimeDerivatives(double[] currentCounts, double simTime){
        double[] rVec = new double[numReactionChannels];

        double[] compartmentDerivatives = getCompartmentDerivatives(currentCounts, simTime);
        // compartment derivatives
        double dSS = compartmentDerivatives[0];
        double dSI = compartmentDerivatives[1];
        double dSC = compartmentDerivatives[2];
        double dSR = compartmentDerivatives[3];
        double dIS = compartmentDerivatives[4];
        double dII = compartmentDerivatives[5];
        double dIC = compartmentDerivatives[6];
        double dIR = compartmentDerivatives[7];
        double dCS = compartmentDerivatives[8];
        double dCI = compartmentDerivatives[9];
        double dCC = compartmentDerivatives[10];
        double dCR = compartmentDerivatives[11];
        double dRS = compartmentDerivatives[12];
        double dRI = compartmentDerivatives[13];
        double dRC = compartmentDerivatives[14];
        double dRR = compartmentDerivatives[15];

        // transmission rate derivatives
        double dBetaOne = getTranRateOneDerivative(simTime);
        double dBetaTwo = getTranRateTwoDerivative(simTime);

        // compartment counts
        double numSS = currentCounts[0];
        double numSI = currentCounts[1];
        double numSC = currentCounts[2];
        double numSR = currentCounts[3];
        double numIS = currentCounts[4];
        double numII = currentCounts[5];
        double numIC = currentCounts[6];
        double numIR = currentCounts[7];
        double numCS = currentCounts[8];
        double numCI = currentCounts[9];
        double numCC = currentCounts[10];
        double numCR = currentCounts[11];
        double numRS = currentCounts[12];
        double numRI = currentCounts[13];
        double numRC = currentCounts[14];
        double numRR = currentCounts[15];

        double transmissionRateOne = getTransmissionRateOne(simTime);
        double moveToCRateOne = rateParameters.get(1).getParameterValue(0);
        double moveToRRateOne = rateParameters.get(2).getParameterValue(0);
        double resusRateOne = rateParameters.get(4).getParameterValue(0);
        double transmissionRateTwo = getTransmissionRateTwo(simTime);
        double moveToCRateTwo = rateParameters.get(6).getParameterValue(0);
        double moveToRRateTwo = rateParameters.get(7).getParameterValue(0);
        double resusRateTwo = rateParameters.get(9).getParameterValue(0);

        double infectionRateModulationI = rateParameters.get(10).getParameterValue(0);
        double infectionRateModulationC = rateParameters.get(11).getParameterValue(0);
        double infectionRateModulationR = rateParameters.get(12).getParameterValue(0);

        double birthMortalityRate = rateParameters.get(13).getParameterValue(0);
        double importationRateOne = rateParameters.get(14).getParameterValue(0);
        double importationRateTwo = rateParameters.get(15).getParameterValue(0);

        // Infection Reactions

        // SS + IS -> 2IS
        rVec[0] = dBetaOne*numSS*numIS + transmissionRateOne*dSS*numIS + transmissionRateOne*numSS*dIS;
        // SS + II -> IS + II
        rVec[1] = dBetaOne*numSS*numII + transmissionRateOne*dSS*numII + transmissionRateOne*numSS*dII;
        // SS + IC -> IS + IC
        rVec[2] = dBetaOne*numSS*numIC + transmissionRateOne*dSS*numIC + transmissionRateOne*numSS*dIC;
        // SS + IR -> IS + IR
        rVec[3] = dBetaOne*numSS*numIR + transmissionRateOne*dSS*numIR + transmissionRateOne*numSS*dIR;

        // SI + IS -> II + IS
        rVec[4] = infectionRateModulationI*(dBetaOne*numSI*numIS + transmissionRateOne*dSI*numIS + transmissionRateOne*numSI*dIS);
        // SI + II -> 2II
        rVec[5] = infectionRateModulationI*(dBetaOne*numSI*numII + transmissionRateOne*dSI*numII + transmissionRateOne*numSI*dII);
        // SI + IC -> II + IC
        rVec[6] = infectionRateModulationI*(dBetaOne*numSI*numIC + transmissionRateOne*dSI*numIC + transmissionRateOne*numSI*dIC);
        // SI + IR -> II + IR
        rVec[7] = infectionRateModulationI*(dBetaOne*numSI*numIR + transmissionRateOne*dSI*numIR + transmissionRateOne*numSI*dIR);

        // SC + IS -> IC + IS
        rVec[8] = infectionRateModulationC*(dBetaOne*numSC*numIS + transmissionRateOne*dSC*numIS + transmissionRateOne*numSC*dIS);
        // SC + II -> IC + II
        rVec[9] = infectionRateModulationC*(dBetaOne*numSC*numII + transmissionRateOne*dSC*numII + transmissionRateOne*numSC*dII);
        // SC + IC -> 2IC
        rVec[10] = infectionRateModulationC*(dBetaOne*numSC*numIC + transmissionRateOne*dSC*numIC + transmissionRateOne*numSC*dIC);
        // SC + IR -> IC + IR
        rVec[11] = infectionRateModulationC*(dBetaOne*numSC*numIR + transmissionRateOne*dSC*numIR + transmissionRateOne*numSC*dIR);

        // SR + IS -> IR + IS
        rVec[12] = infectionRateModulationR*(dBetaOne*numSR*numIS + transmissionRateOne*dSR*numIS + transmissionRateOne*numSR*dIS);
        // SR + II -> IR + II
        rVec[13] = infectionRateModulationR*(dBetaOne*numSR*numII + transmissionRateOne*dSR*numII + transmissionRateOne*numSR*dII);
        // SR + IC -> IR + IC
        rVec[14] = infectionRateModulationR*(dBetaOne*numSR*numIC + transmissionRateOne*dSR*numIC + transmissionRateOne*numSR*dIC);
        // SR + IR -> 2IR
        rVec[15] = infectionRateModulationR*(dBetaOne*numSR*numIR + transmissionRateOne*dSR*numIR + transmissionRateOne*numSR*dIR);

        // SS + SI -> 2SI
        rVec[16] = dBetaTwo*numSS*numSI + transmissionRateTwo*dSS*numSI + transmissionRateTwo*numSS*dSI;
        // SS + II -> SI + II
        rVec[17] = dBetaTwo*numSS*numII + transmissionRateTwo*dSS*numII + transmissionRateTwo*numSS*dII;
        // SS + CI -> SI + CI
        rVec[18] = dBetaTwo*numSS*numCI + transmissionRateTwo*dSS*numCI + transmissionRateTwo*numSS*dCI;
        // SS + RI -> SI + RI
        rVec[19] = dBetaTwo*numSS*numRI + transmissionRateTwo*dSS*numRI + transmissionRateTwo*numSS*dRI;

        // IS + SI -> II + SI
        rVec[20] = infectionRateModulationI*(dBetaTwo*numIS*numSI + transmissionRateTwo*dIS*numSI + transmissionRateTwo*numIS*dSI);
        // IS + II -> 2II
        rVec[21] = infectionRateModulationI*(dBetaTwo*numIS*numII + transmissionRateTwo*dIS*numII + transmissionRateTwo*numIS*dII);
        // IS + CI -> II + CI
        rVec[22] = infectionRateModulationI*(dBetaTwo*numIS*numCI + transmissionRateTwo*dIS*numCI + transmissionRateTwo*numIS*dCI);
        // IS + RI -> II + RI
        rVec[23] = infectionRateModulationI*(dBetaTwo*numIS*numRI + transmissionRateTwo*dIS*numRI + transmissionRateTwo*numIS*dRI);

        // CS + SI -> CI + SI
        rVec[24] = infectionRateModulationC*(dBetaTwo*numCS*numSI + transmissionRateTwo*dCS*numSI + transmissionRateTwo*numCS*dSI);
        // CS + II -> CI + II
        rVec[25] = infectionRateModulationC*(dBetaTwo*numCS*numII + transmissionRateTwo*dCS*numII + transmissionRateTwo*numCS*dII);
        // CS + CI -> 2CI
        rVec[26] = infectionRateModulationC*(dBetaTwo*numCS*numCI + transmissionRateTwo*dCS*numCI + transmissionRateTwo*numCS*dCI);
        // CS + RI -> CI + RI
        rVec[27] = infectionRateModulationC*(dBetaTwo*numCS*numRI + transmissionRateTwo*dCS*numRI + transmissionRateTwo*numCS*dRI);

        // RS + SI -> RI + SI
        rVec[28] = infectionRateModulationR*(dBetaTwo*numRS*numSI + transmissionRateTwo*dRS*numSI + transmissionRateTwo*numRS*dSI);
        // RS + II -> RI + II
        rVec[29] = infectionRateModulationR*(dBetaTwo*numRS*numII + transmissionRateTwo*dRS*numII + transmissionRateTwo*numRS*dII);
        // RS + CI -> RI + CI
        rVec[30] = infectionRateModulationR*(dBetaTwo*numRS*numCI + transmissionRateTwo*dRS*numCI + transmissionRateTwo*numRS*dCI);
        // RS + RI -> 2RI
        rVec[31] = infectionRateModulationR*(dBetaTwo*numRS*numRI + transmissionRateTwo*dRS*numRI + transmissionRateTwo*numRS*dRI);

        // Recovery reactions

        // IS -> CS
        rVec[32] = moveToCRateOne*dIS;
        // II -> CI
        rVec[33] = moveToCRateOne*dII;
        // IC -> CC
        rVec[34] = moveToCRateOne*dIC;
        // IR -> CR
        rVec[35] = moveToCRateOne*dIR;
        // SI -> SC
        rVec[36] = moveToCRateTwo*dSI;
        // II -> IC
        rVec[37] = moveToCRateTwo*dII;
        // CI -> CC
        rVec[38] = moveToCRateTwo*dCI;
        // RI -> RC
        rVec[39] = moveToCRateTwo*dRI;

        // Loss of cross-immunity reactions

        // CS -> RS
        rVec[40] = moveToRRateOne*dCS;
        // CI -> RI
        rVec[41] = moveToRRateOne*dCI;
        // CC -> RC
        rVec[42] = moveToRRateOne*dCC;
        // CR -> RR
        rVec[43] = moveToRRateOne*dCR;
        // SC -> SR
        rVec[44] = moveToRRateTwo*dSC;
        // IC -> IR
        rVec[45] = moveToRRateTwo*dIC;
        // CC -> CR
        rVec[46] = moveToRRateTwo*dCC;
        // RC -> RR
        rVec[47] = moveToRRateTwo*dRC;

        // Resusceptibility reactions

        // RS -> SS
        rVec[48] = resusRateOne*dRS;
        // RI -> SI
        rVec[49] = resusRateOne*dRI;
        // RC -> SC
        rVec[50] = resusRateOne*dRC;
        // RR -> SR
        rVec[51] = resusRateOne*dRR;
        // SR -> SS
        rVec[52] = resusRateTwo*dSR;
        // IR -> IS
        rVec[53] = resusRateTwo*dIR;
        // CR -> CS
        rVec[54] = resusRateTwo*dCR;
        // RR -> RS
        rVec[55] = resusRateTwo*dRR;

        // Infection due to importation

        // SS -> IS
        rVec[56] = dBetaOne*importationRateOne*numSS + transmissionRateOne*importationRateOne*dSS;
        // SI -> II
        rVec[57] = infectionRateModulationI*(dBetaOne*importationRateOne*numSI + transmissionRateOne*importationRateOne*dSI);
        // SC -> IC
        rVec[58] = infectionRateModulationC*(dBetaOne*importationRateOne*numSC + transmissionRateOne*importationRateOne*dSC);
        // SR -> IR
        rVec[59] = infectionRateModulationR*(dBetaOne*importationRateOne*numSR + transmissionRateOne*importationRateOne*dSR);
        // SS -> SI
        rVec[60] = dBetaTwo*importationRateTwo*numSS + transmissionRateTwo*importationRateTwo*dSS;
        // IS -> II
        rVec[61] = infectionRateModulationI*(dBetaTwo*importationRateTwo*numIS + transmissionRateTwo*importationRateTwo*dIS);
        // CS -> CI
        rVec[62] = infectionRateModulationC*(dBetaTwo*importationRateTwo*numCS + transmissionRateTwo*importationRateTwo*dCS);
        // RS -> RI
        rVec[63] = infectionRateModulationR*(dBetaTwo*importationRateTwo*numRS + transmissionRateTwo*importationRateTwo*dRS);


        // Demographic replacement (birth-death reactions)
        // NOTE: we are omitting the trivial SS -> 0, 0 -> SS reaction, so there are 79 reactions in total
        // rather than the 80 that there should be in theory

        // IS -> 0, 0 -> SS
        rVec[64] = birthMortalityRate*dIS;
        // CS -> 0, 0 -> SS
        rVec[65] = birthMortalityRate*dCS;
        // RS -> 0, 0 -> SS
        rVec[66] = birthMortalityRate*dRS;
        // SI -> 0, 0 -> SS
        rVec[67] = birthMortalityRate*dSI;
        // II -> 0, 0 -> SS
        rVec[68] = birthMortalityRate*dII;
        // CI -> 0, 0 -> SS
        rVec[69] = birthMortalityRate*dCI;
        // RI -> 0, 0 -> SS
        rVec[70] = birthMortalityRate*dRI;
        // SC -> 0, 0 -> SS
        rVec[71] = birthMortalityRate*dSC;
        // IC -> 0, 0 -> SS
        rVec[72] = birthMortalityRate*dIC;
        // CC -> 0, 0 -> SS
        rVec[73] = birthMortalityRate*dCC;
        // RC -> 0, 0 -> SS
        rVec[74] = birthMortalityRate*dRC;
        // SR -> 0, 0 -> SS
        rVec[75] = birthMortalityRate*dSR;
        // IR -> 0, 0 -> SS
        rVec[76] = birthMortalityRate*dIR;
        // CR -> 0, 0 -> SS
        rVec[77] = birthMortalityRate*dCR;
        // RR -> 0, 0 -> SS
        rVec[78] = birthMortalityRate*dRR;

        /*
        // SS to IS
        returnVec[0] = transRateOne*dSS*numIS + transRateOne*numSS*dIS;
        returnVec[1] = transRateOne*dSS*numII + transRateOne*numSS*dII;
        returnVec[2] = transRateOne*dSS*numIC + transRateOne*numSS*dIC;
        returnVec[3] = transRateOne*dSS*numIR + transRateOne*numSS*dIR;
        // SS to SI
        returnVec[4] = transRateTwo*dSS*numSI + transRateTwo*numSS*dSI;
        returnVec[5] = transRateTwo*dSS*numII + transRateTwo*numSS*dII;
        returnVec[6] = transRateTwo*dSS*numCI + transRateTwo*numSS*dCI;
        returnVec[7] = transRateTwo*dSS*numRI + transRateTwo*numSS*dRI;
        // SI to II
        returnVec[8] = infectionRateModulationI*transRateOne*dSI*numIS + infectionRateModulationI*transRateOne*numSI*dIS;
        returnVec[9] = infectionRateModulationI*transRateOne*dSI*numII + infectionRateModulationI*transRateOne*numSI*dII;
        returnVec[10] = infectionRateModulationI*transRateOne*dSI*numIC + infectionRateModulationI*transRateOne*numSI*dIC;
        returnVec[11] = infectionRateModulationI*transRateOne*dSI*numIR + infectionRateModulationI*transRateOne*numSI*dIR;
        // SI to SC
        returnVec[12] = moveToCRateTwo*dSI;
        // SC to IC
        returnVec[13] = infectionRateModulationC*transRateOne*dSC*numIS + infectionRateModulationC*transRateOne*numSC*dIS;
        returnVec[14] = infectionRateModulationC*transRateOne*dSC*numII + infectionRateModulationC*transRateOne*numSC*dII;
        returnVec[15] = infectionRateModulationC*transRateOne*dSC*numIC + infectionRateModulationC*transRateOne*numSC*dIC;
        returnVec[16] = infectionRateModulationC*transRateOne*dSC*numIR + infectionRateModulationC*transRateOne*numSC*dIR;
        // SC to SR
        returnVec[17] = moveToRRateTwo*dSC;
        // SR to IR
        returnVec[18] = transRateOne*dSR*numIS + transRateOne*numSR*dIS;
        returnVec[19] = transRateOne*dSR*numII + transRateOne*numSR*dII;
        returnVec[20] = transRateOne*dSR*numIC + transRateOne*numSR*dIC;
        returnVec[21] = transRateOne*dSR*numIR + transRateOne*numSR*dIR;
        // SR to SS
        returnVec[22] = resusRateTwo*dSR;
        // IS to CS
        returnVec[23] = moveToCRateOne*dIS;
        // IS to II
        returnVec[24] = infectionRateModulationI*transRateTwo*dIS*numSI + infectionRateModulationI*transRateTwo*numIS*dSI;
        returnVec[25] = infectionRateModulationI*transRateTwo*dIS*numII + infectionRateModulationI*transRateTwo*numIS*dII;
        returnVec[26] = infectionRateModulationI*transRateTwo*dIS*numCI + infectionRateModulationI*transRateTwo*numIS*dCI;
        returnVec[27] = infectionRateModulationI*transRateTwo*dIS*numRI + infectionRateModulationI*transRateTwo*numIS*dRI;
        // II to CI
        returnVec[28] = infectionRateModulationR*moveToCRateOne*dII;
        // II to IC
        returnVec[29] = infectionRateModulationR*moveToCRateTwo*dII;
        // IC to CC
        returnVec[30] = moveToCRateOne*dIC;
        // IC to IR
        returnVec[31] = moveToRRateTwo*dIC;
        // IR to CR
        returnVec[32] = moveToCRateOne*dIR;
        // IR to IS
        returnVec[33] = resusRateTwo*dIR;
        // CS to RS
        returnVec[34] = moveToRRateOne*dCS;
        // CS to CI
        returnVec[35] = infectionRateModulationC*transRateTwo*dCS*numSI + infectionRateModulationC*transRateTwo*numCS*dSI;
        returnVec[36] = infectionRateModulationC*transRateTwo*dCS*numII + infectionRateModulationC*transRateTwo*numCS*dII;
        returnVec[37] = infectionRateModulationC*transRateTwo*dCS*numCI + infectionRateModulationC*transRateTwo*numCS*dCI;
        returnVec[38] = infectionRateModulationC*transRateTwo*dCS*numRI + infectionRateModulationC*transRateTwo*numCS*dRI;
        // CI to RI
        returnVec[39] = moveToRRateOne*dCI;
        // CI to CC
        returnVec[40] = moveToCRateTwo*dCI;
        // CC to RC
        returnVec[41] = moveToRRateOne*dCC;
        // CC to CR
        returnVec[42] = moveToRRateTwo*dCC;
        // CR to RR
        returnVec[43] = moveToRRateOne*dCR;
        // CR to CS
        returnVec[44] = resusRateTwo*dCR;
        // RS to SS
        returnVec[45] = resusRateOne*dRS;
        // RS to RI
        returnVec[46] = transRateTwo*dRS*numSI + transRateTwo*numRS*dSI;
        returnVec[47] = transRateTwo*dRS*numII + transRateTwo*numRS*dII;
        returnVec[48] = transRateTwo*dRS*numCI + transRateTwo*numRS*dCI;
        returnVec[49] = transRateTwo*dRS*numRI + transRateTwo*numRS*dRI;
        // RI to SI
        returnVec[50] = resusRateOne*dRI;
        // RI to RC
        returnVec[51] = moveToCRateTwo*dRI;
        // RC to SC
        returnVec[52] = resusRateOne*dRC;
        // RC to RR
        returnVec[53] = moveToRRateTwo*dRC;
        // RR to SR
        returnVec[54] = resusRateOne*dRR;
        // RR to RS
        returnVec[55] = resusRateTwo*dRR;
        */
        return rVec;
    }


    protected double[] getSALPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau, double simTime){
        double[] returnVal = new double[numReactionChannels];
        // for standard tau leaping
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = reactionInt[r]*tau;
        }
        // for SAL algorithm, need to add extra terms to account for linear change in intensity
        double[] timeDerivatives = getTimeDerivatives(currentCounts, simTime);
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = returnVal[r] + timeDerivatives[r]*tau*tau*0.5;
        }
        return returnVal;
    }
    // getSALPoissonIntensities should have dimension 79

    // rCounts has reaction counts
    protected double[] getUpdatedCompartmentCounts(double[] currentCounts, double[] rCounts){
        //throw new RuntimeException("getUpdatedCompartmentCounts is definitely being called");

        double[] updatedCounts = new double[numSpecies];

        for (int i = 0; i < numSpecies; i++) {
            updatedCounts[i] = currentCounts[i];
            for (int j = 0; j < numReactionChannels; j++) {
                updatedCounts[i] += vMatrix[i][j]*rCounts[j];
            }
        }

        /*
        // SS(t+tau) = SS(t) - rxn0 - rxn1 - rxn2 - rxn3 - rxn4 - rxn5 - rxn6 - rxn7 + rxn22 + rxn45
        updatedCounts[0] = currentCounts[0] - countsNew[0] - countsNew[1] - countsNew[2] - countsNew[3] - countsNew[4] - countsNew[5] - countsNew[6] - countsNew[7] + countsNew[22] + countsNew[45];
        // SI(t+tau) = SI(t) - rxn8 - rxn9 - rxn10 - rxn11 - rxn12 + rxn4 + rxn5 + rxn6 + rxn7 + rxn50
        updatedCounts[1] = currentCounts[1] - countsNew[8] - countsNew[9] - countsNew[10] - countsNew[11] - countsNew[12] + countsNew[4] + countsNew[5] + countsNew[6] + countsNew[7] + countsNew[50];
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
        */
        return updatedCounts;

        /*
        // ===== DEBUG BLOCK START =====

        System.out.println("----------------------------------------");
        System.out.println("Entering getUpdatedCompartmentCounts()");

        double totalBefore = 0.0;
        double totalAfter = 0.0;

        for (int i = 0; i < currentCounts.length; i++) {
            totalBefore += currentCounts[i];
            totalAfter += updatedCounts[i];
        }

        System.out.println("totalBefore = " + totalBefore);
        System.out.println("totalAfter  = " + totalAfter);

        System.out.println("Current counts:");
        for (int i = 0; i < currentCounts.length; i++) {
            System.out.printf("%2d : %8.1f%n", i, currentCounts[i]);
        }

        System.out.println("Reaction firings:");
        for (int j = 0; j < countsNew.length; j++) {
            if (countsNew[j] != 0.0) {
                System.err.printf("rxn%2d : %8.1f%n", j, countsNew[j]);
            }
        }

        System.out.println("Updated counts:");
        for (int i = 0; i < updatedCounts.length; i++) {
            System.out.printf("%2d : %8.1f%n", i, updatedCounts[i]);
        }

        System.out.println("----------------------------------------");

        // ===== DEBUG BLOCK END =====
        */

        //return updatedCounts;

    }


    protected double getOldestOrigin() {
        if(originOne.getParameterValue(0) >= originTwo.getParameterValue(0)){
            return originOne.getParameterValue(0);
        }else{
            return originTwo.getParameterValue(0);
        }
    }

    @Override
    public double[] introduceSecondPathogen(double simulationTime, double[] currentCounts) {

        // check if second pathogen has not yet been introduced
        if(!secondPathogenIntroduced) {
            double origOne = originOne.getParameterValue(0);
            double origTwo = originTwo.getParameterValue(0);
            // forward time of simulation start time is 0.0, corresponds to backward time of oldest origin
            // in forward time, time of younger origin is origTimeDif
            double origTimeDiff = Math.abs(origOne - origTwo);
            // check if time of younger origin (in forward time) is <= simulationTime
            if(origTimeDiff <= simulationTime){

                currentCounts[0] = currentCounts[0] - 1; // SS

                if (origOne < origTwo) {
                    // pathogen 1 is younger, introduce pathogen 1
                    currentCounts[4] = currentCounts[4] + 1; // IS
                } else {
                    // pathogen 2 is younger, introduce pathogen 2
                    currentCounts[1] = currentCounts[1] + 1; // SI
                }
                secondPathogenIntroduced = true;
                System.out.println("Both pathogens now active");
            }
        }
        return currentCounts;
    }

    protected double getForceOfInfectionOne(double[] currentCounts, double simTime) {
        double tRateOne = getTransmissionRateOne(simTime);
        double numIS = currentCounts[4];
        double numII = currentCounts[5];
        double numIC = currentCounts[6];
        double numIR = currentCounts[7];
        double importationRateOne = rateParameters.get(14).getParameterValue(0);

        return tRateOne*(numIS + numII + numIC + numIR + importationRateOne);
    }

    protected double getForceOfInfectionTwo(double[] currentCounts, double simTime) {
        double tRateTwo = getTransmissionRateTwo(simTime);
        double numSI = currentCounts[1];
        double numII = currentCounts[5];
        double numCI = currentCounts[9];
        double numRI = currentCounts[13];
        double importationRateTwo = rateParameters.get(15).getParameterValue(0);

        return tRateTwo*(numSI + numII + numCI + numRI + importationRateTwo);
    }

    protected double getTransmissionRateOne(double simTime) {
        // NEED TO DO: compute an offset for simTime, if necessary
        double tRateOne = rateParameters.get(0).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            tRateOne = tRateOne*origTimeNumSS*(1.0 + seasonalAmpOne*Math.cos(2.0*Math.PI*(simTime-seasonalPeakDayOne)/seasonalPeriod));
        }
        return tRateOne;
    }

    protected double getTransmissionRateTwo(double simTime) {
        // NEED TO DO: compute an offset for simTime, if necessary
        double tRateTwo = rateParameters.get(5).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            tRateTwo = tRateTwo*origTimeNumSS*(1.0 + seasonalAmpTwo*Math.cos(2.0*Math.PI*(simTime
                    -seasonalPeakDayTwo)/seasonalPeriod));
        }
        return tRateTwo;
    }

    protected double getTranRateOneDerivative(double simTime) {
        // NEED TO DO: compute an offset for simTime, if necessary
        double returnVal = 0.0;
        double tRateOne = rateParameters.get(0).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            returnVal = -tRateOne*origTimeNumSS*seasonalAmpOne*(2.0*Math.PI/seasonalPeriod)
                    *Math.sin(2.0*Math.PI*(simTime-seasonalPeakDayOne)/seasonalPeriod);
        }
        return returnVal;
    }

    protected double getTranRateTwoDerivative(double simTime) {
        // NEED TO DO: compute an offset for simTime, if necessary
        double returnVal = 0.0;
        double tRateTwo = rateParameters.get(5).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            returnVal = -tRateTwo*origTimeNumSS*seasonalAmpTwo*(2.0*Math.PI/seasonalPeriod)
                    *Math.sin(2.0*Math.PI*(simTime-seasonalPeakDayTwo)/seasonalPeriod);
        }
        return returnVal;
    }
}