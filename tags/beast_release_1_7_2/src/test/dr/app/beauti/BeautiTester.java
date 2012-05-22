/*
 * BeautiTester.java
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

package test.dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: BeautiTester.java,v 1.2 2005/07/11 14:07:25 rambaut Exp $
 */
public class BeautiTester {     

    public BeautiTester() {   
        BeautiTesterConfig btc = new BeautiTesterConfig();
        btc.createScriptWriter("tests/run_script.sh");
    	
        BeautiOptions beautiOptions = btc.createOptions();
        btc.importFromFile("examples/Primates.nex", beautiOptions, false); 
        btc.buildNucModels("tests/pri_", beautiOptions);

        beautiOptions = btc.createOptions();
        btc.importFromFile("examples/Primates.nex", beautiOptions, true);
        btc.buildAAModels("tests/pri_", beautiOptions);

        beautiOptions = btc.createOptions();
        btc.importFromFile("examples/Dengue4.env.nex", beautiOptions, false);
//        beautiOptions.fixedSubstitutionRate = false;
        btc.buildNucModels("tests/den_", beautiOptions);

        beautiOptions = btc.createOptions();
        btc.importFromFile("examples/Dengue4.env.nex", beautiOptions, true);
//        beautiOptions.fixedSubstitutionRate = false;
        btc.buildAAModels("tests/den_", beautiOptions);

        btc.closeScriptWriter();
    }

   

    //Main method
    public static void main(String[] args) {

        new BeautiTester();
    }
}

