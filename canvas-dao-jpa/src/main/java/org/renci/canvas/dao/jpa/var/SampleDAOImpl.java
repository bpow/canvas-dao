package org.renci.canvas.dao.jpa.var;

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
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.renci.canvas.dao.var.SampleDAO;
import org.renci.canvas.dao.var.model.Project;
import org.renci.canvas.dao.var.model.Project_;
import org.renci.canvas.dao.var.model.Sample;
import org.renci.canvas.dao.var.model.Sample_;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import org.springframework.transaction.annotation.Transactional;

@Component
@org.springframework.transaction.annotation.Transactional
@OsgiServiceProvider(classes = { SampleDAO.class })
@javax.transaction.Transactional
@Singleton
public class SampleDAOImpl extends BaseDAOImpl<Sample, Integer> implements SampleDAO {

    private static final Logger logger = LoggerFactory.getLogger(SampleDAOImpl.class);

    public SampleDAOImpl() {
        super();
    }

    @Override
    public Class<Sample> getPersistentClass() {
        return Sample.class;
    }

    @org.springframework.transaction.annotation.Transactional
    @javax.transaction.Transactional
    @Override
    public void delete(Sample entity) throws CANVASDAOException {
        logger.debug("ENTERING delete(T)");
        Sample foundEntity = getEntityManager().find(getPersistentClass(), entity.getName());
        getEntityManager().remove(foundEntity);
    }

    @Override
    @Transactional(readOnly = true)
    public List<Sample> findByNameAndProjectName(String name, String projectName) throws CANVASDAOException {
        logger.debug("ENTERING findByNameAndProjectName(String, String)");
        List<Sample> ret = new ArrayList<Sample>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<Sample> crit = critBuilder.createQuery(getPersistentClass());
            Root<Sample> root = crit.from(getPersistentClass());

            List<Predicate> predicates = new ArrayList<Predicate>();
            predicates.add(critBuilder.equal(root.get(Sample_.name), name));
            Join<Sample, Project> sampleProjectJoin = root.join(Sample_.project);
            predicates.add(critBuilder.equal(sampleProjectJoin.get(Project_.id), projectName));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<Sample> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            throw new CANVASDAOException(e);
        }
        return ret;
    }

    @Override
    @Transactional(readOnly = true)
    public List<Sample> findByName(String name) throws CANVASDAOException {
        logger.debug("ENTERING findByName(String)");
        List<Sample> ret = new ArrayList<Sample>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<Sample> crit = critBuilder.createQuery(getPersistentClass());
            Root<Sample> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            predicates.add(critBuilder.equal(root.get(Sample_.name), name));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<Sample> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            throw new CANVASDAOException(e);
        }
        return ret;
    }

}
