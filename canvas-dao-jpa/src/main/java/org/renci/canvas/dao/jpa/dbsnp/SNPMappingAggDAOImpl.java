package org.renci.canvas.dao.jpa.dbsnp;

import java.util.ArrayList;
import java.util.List;

import javax.inject.Singleton;
import javax.persistence.TypedQuery;
import javax.persistence.criteria.CriteriaBuilder;
import javax.persistence.criteria.CriteriaQuery;
import javax.persistence.criteria.Join;
import javax.persistence.criteria.Predicate;
import javax.persistence.criteria.Root;

import org.ops4j.pax.cdi.api.OsgiServiceProvider;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.dbsnp.SNPMappingAggDAO;
import org.renci.canvas.dao.dbsnp.model.SNPMappingAgg;
import org.renci.canvas.dao.dbsnp.model.SNPMappingAggPK;
import org.renci.canvas.dao.dbsnp.model.SNPMappingAgg_;
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.renci.canvas.dao.var.model.LocatedVariant_;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional(readOnly = true)
@OsgiServiceProvider(classes = { SNPMappingAggDAO.class })
@javax.transaction.Transactional(javax.transaction.Transactional.TxType.SUPPORTS)
@Singleton
public class SNPMappingAggDAOImpl extends BaseDAOImpl<SNPMappingAgg, SNPMappingAggPK> implements SNPMappingAggDAO {

    private static final Logger logger = LoggerFactory.getLogger(SNPMappingAggDAOImpl.class);

    public SNPMappingAggDAOImpl() {
        super();
    }

    @Override
    public Class<SNPMappingAgg> getPersistentClass() {
        return SNPMappingAgg.class;
    }

    @Override
    public List<SNPMappingAgg> findByLocatedVariantId(Long locatedVariantId) throws CANVASDAOException {
        logger.debug("ENTERING findByLocatedVariantId(Long)");
        List<SNPMappingAgg> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<SNPMappingAgg> crit = critBuilder.createQuery(getPersistentClass());
            Root<SNPMappingAgg> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();

            Join<SNPMappingAgg, LocatedVariant> snpMappingAggLocatedVariantJoin = root.join(SNPMappingAgg_.locatedVariant);
            predicates.add(critBuilder.equal(snpMappingAggLocatedVariantJoin.get(LocatedVariant_.id), locatedVariantId));

            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            crit.distinct(true);
            TypedQuery<SNPMappingAgg> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ret;
    }

}
