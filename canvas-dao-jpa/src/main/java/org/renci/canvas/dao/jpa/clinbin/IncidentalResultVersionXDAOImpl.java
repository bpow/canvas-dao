package org.renci.canvas.dao.jpa.clinbin;

import javax.inject.Singleton;
import javax.persistence.TypedQuery;
import javax.persistence.criteria.CriteriaBuilder;
import javax.persistence.criteria.CriteriaQuery;
import javax.persistence.criteria.Root;

import org.ops4j.pax.cdi.api.OsgiServiceProvider;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.IncidentalResultVersionXDAO;
import org.renci.canvas.dao.clinbin.model.IncidentalResultVersionX;
import org.renci.canvas.dao.clinbin.model.IncidentalResultVersionX_;
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional(readOnly = true)
@OsgiServiceProvider(classes = { IncidentalResultVersionXDAO.class })
@javax.transaction.Transactional(javax.transaction.Transactional.TxType.SUPPORTS)
@Singleton
public class IncidentalResultVersionXDAOImpl extends BaseDAOImpl<IncidentalResultVersionX, Integer> implements IncidentalResultVersionXDAO {

    private static final Logger logger = LoggerFactory.getLogger(IncidentalResultVersionXDAOImpl.class);

    public IncidentalResultVersionXDAOImpl() {
        super();
    }

    @Override
    public Class<IncidentalResultVersionX> getPersistentClass() {
        return IncidentalResultVersionX.class;
    }

    @Override
    public Integer findMaxResultVersion() throws CANVASDAOException {
        logger.debug("ENTERING findMaxResultVersion()");
        CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
        CriteriaQuery<Integer> crit = critBuilder.createQuery(Integer.class);
        Root<IncidentalResultVersionX> root = crit.from(getPersistentClass());
        crit.select(critBuilder.max(root.get(IncidentalResultVersionX_.id)));
        TypedQuery<Integer> query = getEntityManager().createQuery(crit);
        Integer ret = query.getSingleResult();
        return ret;
    }

}
