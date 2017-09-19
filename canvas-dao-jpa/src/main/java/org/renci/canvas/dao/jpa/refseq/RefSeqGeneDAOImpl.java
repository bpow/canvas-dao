package org.renci.canvas.dao.jpa.refseq;

import java.util.ArrayList;
import java.util.List;

import javax.inject.Singleton;
import javax.persistence.TypedQuery;
import javax.persistence.criteria.CriteriaBuilder;
import javax.persistence.criteria.CriteriaQuery;
import javax.persistence.criteria.Join;
import javax.persistence.criteria.Predicate;
import javax.persistence.criteria.Root;

import org.apache.commons.lang3.StringUtils;
import org.ops4j.pax.cdi.api.OsgiServiceProvider;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.renci.canvas.dao.refseq.RefSeqGeneDAO;
import org.renci.canvas.dao.refseq.model.RefSeqGene;
import org.renci.canvas.dao.refseq.model.RefSeqGene_;
import org.renci.canvas.dao.refseq.model.RegionGroup;
import org.renci.canvas.dao.refseq.model.RegionGroup_;
import org.renci.canvas.dao.refseq.model.Transcript;
import org.renci.canvas.dao.refseq.model.Transcript_;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional(readOnly = true)
@OsgiServiceProvider(classes = { RefSeqGeneDAO.class })
@javax.transaction.Transactional(javax.transaction.Transactional.TxType.SUPPORTS)
@Singleton
public class RefSeqGeneDAOImpl extends BaseDAOImpl<RefSeqGene, Integer> implements RefSeqGeneDAO {

    private static final Logger logger = LoggerFactory.getLogger(RefSeqGeneDAOImpl.class);

    public RefSeqGeneDAOImpl() {
        super();
    }

    @Override
    public Class<RefSeqGene> getPersistentClass() {
        return RefSeqGene.class;
    }

    @Override
    public List<RefSeqGene> findByExample(RefSeqGene refseqGene) throws CANVASDAOException {
        logger.debug("ENTERING findByExample(RefSeqGene)");
        List<RefSeqGene> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<RefSeqGene> crit = critBuilder.createQuery(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Root<RefSeqGene> fromRefSeqGene = crit.from(getPersistentClass());

            if (StringUtils.isNotEmpty(refseqGene.getName())) {
                predicates.add(critBuilder.equal(fromRefSeqGene.get(RefSeqGene_.name), refseqGene.getName()));
            }

            if (StringUtils.isNotEmpty(refseqGene.getRefseqVersion())) {
                predicates.add(critBuilder.equal(fromRefSeqGene.get(RefSeqGene_.refseqVersion), refseqGene.getRefseqVersion()));
            }

            if (StringUtils.isNotEmpty(refseqGene.getDescription())) {
                predicates.add(critBuilder.equal(fromRefSeqGene.get(RefSeqGene_.description), refseqGene.getDescription()));
            }
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<RefSeqGene> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ret;
    }

    @Override
    public List<RefSeqGene> findByRefSeqVersion(String refSeqVersion) throws CANVASDAOException {
        logger.debug("ENTERING findByRefSeqVersion(String)");
        List<RefSeqGene> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<RefSeqGene> crit = critBuilder.createQuery(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Root<RefSeqGene> fromRefSeqGene = crit.from(getPersistentClass());
            predicates.add(critBuilder.equal(fromRefSeqGene.get(RefSeqGene_.refseqVersion), refSeqVersion));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            crit.orderBy(critBuilder.asc(fromRefSeqGene.get(RefSeqGene_.name)));
            crit.distinct(true);
            TypedQuery<RefSeqGene> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ret;
    }

    @Override
    public List<RefSeqGene> findByRefSeqVersionAndTranscriptId(String refSeqVersion, String transcriptId) throws CANVASDAOException {
        logger.debug("ENTERING findByRefSeqVersionAndTranscriptId(String, String)");
        List<RefSeqGene> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<RefSeqGene> crit = critBuilder.createQuery(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Root<RefSeqGene> fromRefSeqGene = crit.from(getPersistentClass());

            predicates.add(critBuilder.equal(fromRefSeqGene.get(RefSeqGene_.refseqVersion), refSeqVersion));

            Join<RefSeqGene, RegionGroup> refSeqGeneRegionGroupJoin = fromRefSeqGene.join(RefSeqGene_.locations);
            Join<RegionGroup, Transcript> regionGroupTranscriptJoin = refSeqGeneRegionGroupJoin.join(RegionGroup_.transcript);
            predicates.add(critBuilder.equal(regionGroupTranscriptJoin.get(Transcript_.id), transcriptId));

            crit.orderBy(critBuilder.asc(fromRefSeqGene.get(RefSeqGene_.name)));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            crit.distinct(true);
            TypedQuery<RefSeqGene> query = getEntityManager().createQuery(crit);
            // query.setHint(QueryHints.FETCHGRAPH, getEntityManager().getEntityGraph("withLocations"));
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ret;
    }

    @Override
    public List<RefSeqGene> findByTranscriptId(String transcriptId) throws CANVASDAOException {
        logger.debug("ENTERING findByTranscriptId(String)");
        List<RefSeqGene> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<RefSeqGene> crit = critBuilder.createQuery(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Root<RefSeqGene> fromRefSeqGene = crit.from(getPersistentClass());

            Join<RefSeqGene, RegionGroup> refSeqGeneRegionGroupJoin = fromRefSeqGene.join(RefSeqGene_.locations);
            Join<RegionGroup, Transcript> regionGroupTranscriptJoin = refSeqGeneRegionGroupJoin.join(RegionGroup_.transcript);
            predicates.add(critBuilder.equal(regionGroupTranscriptJoin.get(Transcript_.id), transcriptId));

            crit.orderBy(critBuilder.asc(fromRefSeqGene.get(RefSeqGene_.name)));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            crit.distinct(true);
            TypedQuery<RefSeqGene> query = getEntityManager().createQuery(crit);
            // query.setHint(QueryHints.FETCHGRAPH, getEntityManager().getEntityGraph("withLocations"));
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            e.printStackTrace();
        }
        return ret;
    }

}
