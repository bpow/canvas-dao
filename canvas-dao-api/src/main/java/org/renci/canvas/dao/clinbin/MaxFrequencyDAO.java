package org.renci.canvas.dao.clinbin;

import java.util.List;

import org.renci.canvas.dao.BaseDAO;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.model.MaxFrequency;
import org.renci.canvas.dao.clinbin.model.MaxFrequencyPK;

public interface MaxFrequencyDAO extends BaseDAO<MaxFrequency, MaxFrequencyPK> {

    public List<MaxFrequency> findByGeneNameAndMaxAlleleFrequency(String name, Double threshold) throws CANVASDAOException;

    public List<MaxFrequency> findByLocatedVariantId(Long locatedVariantId) throws CANVASDAOException;

}
