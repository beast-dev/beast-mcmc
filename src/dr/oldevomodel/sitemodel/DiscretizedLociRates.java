/*
 * DiscretizedLociRates.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodel.sitemodel;

import dr.inference.model.*;
import dr.inference.distribution.ParametricDistributionModel;

/**
 * @author Chieh-Hsi Wu
 *
 * This class models the relative loci rates using the given discritized parametric distribution. 
 */
public class DiscretizedLociRates extends AbstractModel {

    private CompoundParameter lociRates;
    private Parameter rateCategoryParameter;
    private ParametricDistributionModel distrModel;
    private double normalizeRateTo;
    private double[] rates;
    private boolean normalize;
    private int categoryCount;
    private double scaleFactor;
    private boolean completeSetup;


    public DiscretizedLociRates(
            CompoundParameter lociRates,
            Parameter rateCategoryParameter,
            ParametricDistributionModel model,
            boolean normalize,
            double normalizeLociRateTo,
            int categoryCount) {
        super("DiscretizedLociRatesModel");
        this.lociRates = lociRates;

        this.rateCategoryParameter = rateCategoryParameter;
        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        this.rateCategoryParameter.addBounds(bound);
        
        this.distrModel = model;
        this.normalizeRateTo = normalizeLociRateTo;
        this.normalize = normalize;
        this.categoryCount = categoryCount;
        rates = new double[categoryCount];
        completeSetup = true;
        setupRates();


        addModel(distrModel);
        addVariable(this.rateCategoryParameter);

    }

    private void setupRates(){
        if(completeSetup){

            double categoryIntervalSize = 1.0/categoryCount;
            for(int i = 0; i < categoryCount; i++){
                rates[i]= distrModel.quantile((i+0.5)*categoryIntervalSize);
            }
        }
        if(normalize){
           computeFactor();
        }

        completeSetup = false;
        int lociCount = rateCategoryParameter.getDimension();
        for(int i = 0; i < lociCount; i ++){
            lociRates.setParameterValue(i,rates[(int)rateCategoryParameter.getParameterValue(i)]*scaleFactor);
        }

    }
    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distrModel) {
            completeSetup = true;
            setupRates();
            //System.out.println("speed investigation 1");
            fireModelChanged();
        }else if (model == rateCategoryParameter) {
            //System.out.println("speed investigation 2");
            setupRates();
            fireModelChanged(null, index);
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        //System.out.println("speed investigation 3");
        setupRates();
        fireModelChanged(null, index);
    }

    protected void storeState() {
    }
    protected void acceptState() {
    }
    protected void restoreState() {
        //setupRates();
    }
    private void computeFactor(){
        double sumRates = 0.0;
        int lociCount = rateCategoryParameter.getDimension();
        for(int i = 0; i < lociCount; i++){
            sumRates += rates[(int)rateCategoryParameter.getParameterValue(i)];
        }
        scaleFactor = normalizeRateTo/(sumRates/lociCount);
    }

}
