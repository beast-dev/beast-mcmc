package dr.app.mapper.application;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface MapperFileMenuHandler {

    Action getImportStrainsAction();

    Action getImportMeasurementsAction();

    Action getImportLocationsAction();

    Action getImportTreesAction();

    Action getExportDataAction();

	Action getExportPDFAction();

}