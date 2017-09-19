package org.renci.canvas.dao.jpa.clinbin;

import java.util.ArrayList;
import java.util.List;

import javax.inject.Singleton;
import javax.persistence.TypedQuery;
import javax.persistence.criteria.CriteriaBuilder;
import javax.persistence.criteria.CriteriaQuery;
import javax.persistence.criteria.Predicate;
import javax.persistence.criteria.Root;

import org.ops4j.pax.cdi.api.OsgiServiceProvider;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.MissingHaplotypeDAO;
import org.renci.canvas.dao.clinbin.model.MissingHaplotype;
import org.renci.canvas.dao.clinbin.model.MissingHaplotypePK;
import org.renci.canvas.dao.clinbin.model.MissingHaplotypePK_;
import org.renci.canvas.dao.clinbin.model.MissingHaplotype_;
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional(readOnly = true)
@OsgiServiceProvider(classes = { MissingHaplotypeDAO.class })
@javax.transaction.Transactional(javax.transaction.Transactional.TxType.SUPPORTS)
@Singleton
public class MissingHaplotypeDAOImpl extends BaseDAOImpl<MissingHaplotype, MissingHaplotypePK> implements MissingHaplotypeDAO {

    private static final Logger logger = LoggerFactory.getLogger(MissingHaplotypeDAOImpl.class);

    public MissingHaplotypeDAOImpl() {
        super();
    }

    @Override
    public Class<MissingHaplotype> getPersistentClass() {
        return MissingHaplotype.class;
    }

    @Override
    public List<MissingHaplotype> findByParticipantAndIncidentalBinIdAndListVersion(String participantId, Integer incidentalBinId,
            Integer listVersion) throws CANVASDAOException {
        logger.debug("ENTERING findByParticipantAndIncidentalBinIdAndListVersion(String, Integer, Integer)");
        List<MissingHaplotype> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<MissingHaplotype> crit = critBuilder.createQuery(getPersistentClass());
            Root<MissingHaplotype> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            predicates.add(critBuilder.equal(root.get(MissingHaplotype_.id).get(MissingHaplotypePK_.participant), participantId));
            predicates.add(critBuilder.equal(root.get(MissingHaplotype_.id).get(MissingHaplotypePK_.incidentalBin), incidentalBinId));
            predicates.add(critBuilder.equal(root.get(MissingHaplotype_.id).get(MissingHaplotypePK_.listVersion), listVersion));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            crit.distinct(true);
            TypedQuery<MissingHaplotype> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ret;
    }

}
