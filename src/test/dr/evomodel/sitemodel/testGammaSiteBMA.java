/*
 * testGammaSiteBMA.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.evomodel.sitemodel;

import junit.framework.TestCase;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.oldevomodel.sitemodel.GammaSiteBMA;
import dr.evolution.datatype.Nucleotides;

/**
 * @author Chieh-Hsi Wu
 *
 */
public class testGammaSiteBMA extends TestCase {
    interface Instance {
        SubstitutionModel getSubstModel ();

        double getMu();

        double getLogitInvar();

        double getLogShape();

        int getCategoryCount();

        Variable<Integer> getModelChoose();

        double[] getCategoryRates();

        double[] getCategoryProportions();
    }

    //Neither gamma shape nor site invariant parameters is included.
    Instance test0 = new Instance(){
        public SubstitutionModel getSubstModel (){

            //Create a JC model
            Parameter kappa = new Parameter.Default(1, 1);
            double[] pi = new double[]{0.25, 0.25, 0.25, 0.25};

            Parameter freqs = new Parameter.Default(pi);
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
            HKY jc = new HKY(kappa, f);
            return jc;
        }

        public double getMu(){
            return 1.5;
        }


        public double getLogitInvar(){
            //The invariant site proportion, p = 0.2.
            return -1.38629436112;
        }

        public double getLogShape(){
            //Gamma shape parameter, alpha = 2.
            return 0.69314718056;
        }

        public int getCategoryCount(){
            return 4;
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{0, 0});
        }

        public double[] getCategoryRates(){
            return new double[]{0.0, 1.5, 1.5, 1.5, 1.5};
        }

        public double[] getCategoryProportions(){
            return new double[]{0, 0.25, 0.25, 0.25, 0.25};
        }

    };

    //Neither gamma shape nor site invariant parameters is included.
    Instance test1 = new Instance(){
        public SubstitutionModel getSubstModel (){

            //Create a JC model
            Parameter kappa = new Parameter.Default(1, 1);
            double[] pi = new double[]{0.25, 0.25, 0.25, 0.25};

            Parameter freqs = new Parameter.Default(pi);
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
            HKY jc = new HKY(kappa, f);
            return jc;
        }

        public double getMu(){
            return 1.5;
        }


        public double getLogitInvar(){
            //The invariant site proportion, p = 0.2.
            return -1.38629436112;
        }

        public double getLogShape(){
            //Gamma shape parameter, alpha = 2.
            return 0.69314718056;
        }

        public int getCategoryCount(){
            return 4;
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{0, 1});
        }

        public double[] getCategoryRates(){
            return new double[]{0.0, 1.875, 1.875, 1.875, 1.875};
        }

        public double[] getCategoryProportions(){
            return new double[]{0.2, 0.2, 0.2, 0.2, 0.2};
        }


    };

    //Neither gamma shape nor site invariant parameters is included.
    Instance test2 = new Instance(){
        public SubstitutionModel getSubstModel (){

            //Create a JC model
            Parameter kappa = new Parameter.Default(1, 1);
            double[] pi = new double[]{0.25, 0.25, 0.25, 0.25};

            Parameter freqs = new Parameter.Default(pi);
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
            HKY jc = new HKY(kappa, f);
            return jc;
        }

        public double getMu(){
            return 1.5;
        }


        public double getLogitInvar(){
            //The invariant site proportion, p = 0.2.
            return -1.38629436112;
        }

        public double getLogShape(){
            //Gamma shape parameter, alpha = 2.
            return 0.69314718056;
        }

        public int getCategoryCount(){
            return 4;
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{1, 1});
        }

        public double[] getCategoryRates(){
            return new double[]{0.0, 0.59824693778, 1.28130225333, 2.07933195980, 3.54111884909};
        }

        public double[] getCategoryProportions(){
            return new double[]{0.2, 0.2, 0.2, 0.2, 0.2};
        }
    };

    

    Instance test3 = new Instance(){
        public SubstitutionModel getSubstModel (){

            //Create a JC model
            Parameter kappa = new Parameter.Default(1, 1);
            double[] pi = new double[]{0.25, 0.25, 0.25, 0.25};

            Parameter freqs = new Parameter.Default(pi);
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
            HKY jc = new HKY(kappa, f);
            return jc;
        }

        public double getMu(){
            return 1.5;
        }


        public double getLogitInvar(){
            //The invariant site proportion, p = 0.2.
            return -1.38629436112;
        }

        public double getLogShape(){
            //Gamma shape parameter, alpha = 2.
            return 0.69314718056;
        }

        public int getCategoryCount(){
            return 8;
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{1, 1});
        }

        public double[] getCategoryRates(){
            return new double[]{0.0,
                    0.387224876120, 0.757770190732, 1.085824675833, 1.426031312891,
                    1.810616868095, 2.285251426765, 2.955667305794, 4.291613343768};
        }

        public double[] getCategoryProportions(){
            return new double[]{0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        }
    };
    Instance[] all = {test0,test1,test2,test3};
    public void testGammaSiteBMA(){

        for(Instance test: all){
            SubstitutionModel substModel = test.getSubstModel();
            Parameter mu = new Parameter.Default(test.getMu());
            Parameter logitInvar = new Parameter.Default(test.getLogitInvar());
            Parameter logShape = new Parameter.Default(test.getLogShape());
            int catCount = test.getCategoryCount();
            Variable<Integer> modelChoose = test.getModelChoose();
            GammaSiteBMA gammaSiteBMA = new GammaSiteBMA(
                    substModel,
                    mu,
                    logitInvar,
                    logShape,
                    catCount,
                    modelChoose
            );

            double[] catRates = gammaSiteBMA.getCategoryRates();
            double[] expectedCatRates = test.getCategoryRates();

            for(int i = 0; i < catRates.length; i++){
                assertEquals(catRates[i], expectedCatRates[i], 8e-10);
            }

            double[] catProps = gammaSiteBMA.getCategoryProportions();
            double[] expectedCatProps = test.getCategoryProportions();
            for(int i = 0; i < catProps.length; i++){
                assertEquals(catProps[i], expectedCatProps[i], 1e-10);
            }

        }


    }


}
