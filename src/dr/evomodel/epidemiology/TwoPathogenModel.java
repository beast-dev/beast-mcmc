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
    protected double mostRecentSamplingDateOne;
    protected double seasonalAmpTwo;
    protected double seasonalPeakDayOne;
    protected double seasonalPeakDayTwo;
    protected double mostRecentSamplingDateTwo;
    protected double seasonalOffset;
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
            Parameter mostRecentSamplingDateOne,
            Parameter seasonalAmpTwo,
            Parameter seasonalPeakDayTwo,
            Parameter mostRecentSamplingDateTwo) {
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
        //this.totalPopSize = originTimeNumSS.getParameterValue(0) + 1;
        this.numGridPoints = numGridPoints;
        this.cutOff = cutOff;
        this.numReactionChannels = numReactionChannels;
        this.numSpecies = compartmentCounts.size();
        this.vMatrix = getVMatrix();
        this.mostRecentSamplingDateOne = mostRecentSamplingDateOne.getParameterValue(0);
        this.mostRecentSamplingDateTwo = mostRecentSamplingDateTwo.getParameterValue(0);

        setTotalPopSize();

        this.seasonalModel = seasonalModel;
        if(seasonalModel) {
            this.seasonalPeriod = seasonalPeriod.getParameterValue(0);
            this.seasonalAmpOne = seasonalAmpOne.getParameterValue(0);
            this.seasonalAmpTwo = seasonalAmpTwo.getParameterValue(0);
            this.seasonalPeakDayOne = seasonalPeakDayOne.getParameterValue(0);
            this.seasonalPeakDayTwo = seasonalPeakDayTwo.getParameterValue(0);
            setSeasonalOffset();
        }
    }

    protected void setTotalPopSize() {
        //double origOne = originOne.getParameterValue(0);
        //double origTwo = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);
        double originTimeSI = originTimeNumSI;
        double originTimeIS = originTimeNumIS;

        double mostRecentSamplingDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
        double forwardOrigOne = cutOff - originOne.getParameterValue(0) -
                (mostRecentSamplingDate - mostRecentSamplingDateOne);
        double forwardOrigTwo = cutOff - originTwo.getParameterValue(0) -
                (mostRecentSamplingDate - mostRecentSamplingDateTwo);

        if (forwardOrigOne < forwardOrigTwo) {
            // total compartment counts at origin time should be originTimeSS + originTimeIS
            totalPopSize = originTimeSS + originTimeIS;
        } else if (forwardOrigTwo < forwardOrigOne) {
            totalPopSize = originTimeSS + originTimeSI;
        } else {
            // total compartment counts at origin time should be originTimeSS + originTimeIS + originTimeSI
            totalPopSize = originTimeSS + originTimeIS + originTimeSI;
        }
    }

    protected void setSeasonalOffset(){
        double mostRecentSamplingDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
        double simStartDate = mostRecentSamplingDate - cutOff;
        seasonalOffset = simStartDate - Math.floor(simStartDate);
    }

    protected void setOriginTimeCompartmentCounts(int index){
        secondPathogenIntroduced = false;

        double origOne = originOne.getParameterValue(0);
        double origTwo = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);
        double originTimeSI = originTimeNumSI;
        double originTimeIS = originTimeNumIS;
        double mostRecentSamplingDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
        // units of time beyond simulation start time (which corresponds to cutOff time) of pathogen one origin
        double forwardOrigOne = cutOff - origOne - (mostRecentSamplingDate - mostRecentSamplingDateOne);
        // units of time beyond simulation start time (which corresponds to cutOff time) of pathogen two origin
        double forwardOrigTwo = cutOff - origTwo - (mostRecentSamplingDate - mostRecentSamplingDateTwo);

        if (forwardOrigOne < 0) {
            throw new RuntimeException("Origin time of pathogen 1 is further back in time than cutOff. " +
                    "Need to increase value of cutOff");
        }
        if (forwardOrigTwo < 0) {
            throw new RuntimeException("Origin time of pathogen 2 is further back in time than cutOff. " +
                    "Need to increase value of cutOff.");
        }

        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 0);
        }

        // SS = originTimeSS
        compartmentCounts.get(0).setParameterValue(index, originTimeSS);

        if (forwardOrigOne < forwardOrigTwo) {
            // pathogen 1 is older, start in IS
            compartmentCounts.get(4).setParameterValue(index, originTimeIS);
            // total compartment counts at origin time should be originTimeSS + originTimeIS
            if(totalPopSize != originTimeSS + originTimeIS){
                throw new RuntimeException("Total pop size mismatch");
            }
        } else if (forwardOrigTwo < forwardOrigOne) {
            // pathogen 2 is older, start in SI
            compartmentCounts.get(1).setParameterValue(index, originTimeSI);
            // total compartment counts at origin time should be originTimeSS + originTimeSI
            if(totalPopSize != originTimeSS + originTimeSI){
                throw new RuntimeException("Total pop size mismatch");
            }
        } else {
            // no need to "introduce" second pathogen while doing forward time simulation
            secondPathogenIntroduced = true;
            // origins equal
            compartmentCounts.get(1).setParameterValue(index, originTimeSI); // SI
            compartmentCounts.get(4).setParameterValue(index, originTimeIS); // IS
            // total compartment counts at origin time should be originTimeSS + originTimeIS + originTimeSI
            if(totalPopSize != originTimeSS + originTimeIS + originTimeSI){
                throw new RuntimeException("Total pop size mismatch");
            }
        }
    }

    protected void setDefaultCompartmentCounts(int index){
        //double origOne = originOne.getParameterValue(0);
        //double origTwo = originTwo.getParameterValue(0);
        double originTimeSS = originTimeNumSS.getParameterValue(0);
        double originTimeSI = originTimeNumSI;
        double originTimeIS = originTimeNumIS;

        double mostRecentSamplingDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
        double forwardOrigOne = cutOff - originOne.getParameterValue(0) -
                (mostRecentSamplingDate - mostRecentSamplingDateOne);
        double forwardOrigTwo = cutOff - originTwo.getParameterValue(0) -
                (mostRecentSamplingDate - mostRecentSamplingDateTwo);

        // initialize everything to 0
        for (int i = 0; i < compartmentCounts.size(); i++) {
            compartmentCounts.get(i).setParameterValue(index, 0);
        }

        if (forwardOrigOne < forwardOrigTwo) {
            // pathogen 1 is older, start in IS
            // total compartment counts at origin time should be originTimeSS + originTimeIS
            // default SS value should be same
            compartmentCounts.get(0).setParameterValue(index, originTimeSS + originTimeIS);
        } else if (forwardOrigTwo < forwardOrigOne) {
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
    // columns are reactions, rows are compartment counts, there are 56 reactions
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

        return updatedCounts;
    }

    // gets time into the past from most recent sampling date to oldest epidemic origin
    protected double getOldestOrigin() {
        double originOneCalTime = mostRecentSamplingDateOne -
                originOne.getParameterValue(0);
        double originTwoCalTime = mostRecentSamplingDateTwo -
                originTwo.getParameterValue(0);
        double oldestOriginCalTime = Math.min(originOneCalTime, originTwoCalTime);
        double moreRecentSamplingDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
        return moreRecentSamplingDate - oldestOriginCalTime;
    }

    @Override
    public double[] introduceSecondPathogen(double simulationTime, double[] currentCounts) {

        // check if second pathogen has not yet been introduced
        if(!secondPathogenIntroduced) {
            double mostRecentSamplingDate = Math.max(mostRecentSamplingDateOne, mostRecentSamplingDateTwo);
            double forwardOrigOne = cutOff - originOne.getParameterValue(0) -
                    (mostRecentSamplingDate - mostRecentSamplingDateOne);
            double forwardOrigTwo = cutOff - originTwo.getParameterValue(0) -
                    (mostRecentSamplingDate - mostRecentSamplingDateTwo);

            // forward time of simulation start time is 0.0,
            double youngerForwardOrigTime = Math.max(forwardOrigOne, forwardOrigTwo);

            // check if time of younger origin (in forward time) is <= simulationTime
            if(youngerForwardOrigTime <= simulationTime){
                // SS decreases by 1
                currentCounts[0] = currentCounts[0] - 1;

                if (forwardOrigOne > forwardOrigTwo) {
                    // pathogen 1 is younger, introduce pathogen 1
                    currentCounts[4] = currentCounts[4] + 1; // IS
                } else {
                    // pathogen 2 is younger, introduce pathogen 2
                    currentCounts[1] = currentCounts[1] + 1; // SI
                }
                secondPathogenIntroduced = true;
                //System.out.println("Both pathogens now active");
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
        double tRateOne = rateParameters.get(0).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            tRateOne = tRateOne*origTimeNumSS*(1.0 + seasonalAmpOne*Math.cos(2.0*Math.PI*(simTime + seasonalOffset
                    - seasonalPeakDayOne)/seasonalPeriod));
        }
        return tRateOne;
    }

    protected double getTransmissionRateTwo(double simTime) {
        double tRateTwo = rateParameters.get(5).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            tRateTwo = tRateTwo*origTimeNumSS*(1.0 + seasonalAmpTwo*Math.cos(2.0*Math.PI*(simTime + seasonalOffset
                    -seasonalPeakDayTwo)/seasonalPeriod));
        }
        return tRateTwo;
    }

    protected double getTranRateOneDerivative(double simTime) {
        double returnVal = 0.0;
        double tRateOne = rateParameters.get(0).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            returnVal = -tRateOne*origTimeNumSS*seasonalAmpOne*(2.0*Math.PI/seasonalPeriod)
                    *Math.sin(2.0*Math.PI*(simTime + seasonalOffset -seasonalPeakDayOne)/seasonalPeriod);
        }
        return returnVal;
    }

    protected double getTranRateTwoDerivative(double simTime) {
        double returnVal = 0.0;
        double tRateTwo = rateParameters.get(5).getParameterValue(0);
        double origTimeNumSS = originTimeNumSS.getParameterValue(0);
        if(seasonalModel){
            returnVal = -tRateTwo*origTimeNumSS*seasonalAmpTwo*(2.0*Math.PI/seasonalPeriod)
                    *Math.sin(2.0*Math.PI*(simTime + seasonalOffset -seasonalPeakDayTwo)/seasonalPeriod);
        }
        return returnVal;
    }
}