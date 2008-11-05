/*
 * HKY.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.newsubstmodel;

import dr.app.beauti.options.NucModelType;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;


/**
 * Hasegawa-Kishino-Yano model of nucleotide evolution
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: HKY.java,v 1.42 2005/09/23 13:17:59 rambaut Exp $
 */
public class HKY extends BaseSubstitutionModel {

    private Parameter kappaParameter = null;

    /**
     * A constructor which allows a more programmatic approach with
     * fixed kappa.
     * @param kappa
     * @param freqModel
     */
    public HKY(double kappa, FrequencyModel freqModel) {
        this(new Parameter.Default(kappa), freqModel);
    }
    /**
     * Constructor
     * @param kappaParameter
     * @param freqModel
     */
    public HKY(Parameter kappaParameter, FrequencyModel freqModel) {

        super("HKY", Nucleotides.INSTANCE, freqModel);

        this.kappaParameter = kappaParameter;
        addParameter(kappaParameter);
        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        addStatistic(tsTvStatistic);
    }

    /**
     * set kappa
     * @param kappa
     */
    public void setKappa(double kappa) {
        kappaParameter.setParameterValue(0, kappa);
        updateMatrix = true;
    }

    /**
     * @return kappa
     */
    public final double getKappa() {
        return kappaParameter.getParameterValue(0);
    }

    /**
     * set ts/tv
     * @param tsTv
     */
    public void setTsTv(double tsTv) {
        double freqA = freqModel.getFrequency(0);
        double freqC = freqModel.getFrequency(1);
        double freqG = freqModel.getFrequency(2);
        double freqT = freqModel.getFrequency(3);
        double freqR = freqA + freqG;
        double freqY = freqC + freqT;
        setKappa((tsTv * freqR * freqY) / (freqA * freqG + freqC * freqT));
    }

    /**
     * @return tsTv
     */
    public double getTsTv() {
        double freqA = freqModel.getFrequency(0);
        double freqC = freqModel.getFrequency(1);
        double freqG = freqModel.getFrequency(2);
        double freqT = freqModel.getFrequency(3);
        double freqR = freqA + freqG;
        double freqY = freqC + freqT;
        double tsTv = (getKappa() * (freqA * freqG + freqC * freqT)) / (freqR * freqY);

        return tsTv;
    }

    protected void frequenciesChanged() {
    }

    protected void ratesChanged() {
    }

    protected void setupRelativeRates(double[] rates) {
        double kappa =  kappaParameter.getParameterValue(0);
        rates[0] = 1.0;
        rates[1] = kappa;
        rates[2] = 1.0;
        rates[3] = 1.0;
        rates[4] = kappa;
        rates[5] = 1.0;
    }

    //
    // Private stuff
    //

    private Statistic tsTvStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "tsTv";
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int dim) {
            return getTsTv();
        }

    };
}