package org.renci.canvas.dao.clinbin;

import java.util.List;

import org.renci.canvas.dao.BaseDAO;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.model.IncidentalBinGeneX;
import org.renci.canvas.dao.clinbin.model.IncidentalBinGeneXPK;

public interface IncidentalBinGeneXDAO extends BaseDAO<IncidentalBinGeneX, IncidentalBinGeneXPK> {

    public List<IncidentalBinGeneX> findByIncidentalBinIdAndVersionAndZygosityModes(Integer incidentalBinId, Integer incidentalBinVersion,
            List<String> zygosityModeList) throws CANVASDAOException;

    public List<IncidentalBinGeneX> findByIncidentalBinIdAndVersion(Integer incidentalBinId, Integer incidentalBinVersion)
            throws CANVASDAOException;

}
