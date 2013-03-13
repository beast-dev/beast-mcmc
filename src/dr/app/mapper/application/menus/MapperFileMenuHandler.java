package dr.app.mapper.application.menus;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface MapperFileMenuHandler {

    Action getImportMeasurementsAction();

    Action getImportLocationsAction();

    Action getImportTreesAction();

    Action getExportDataAction();

	Action getExportPDFAction();

}