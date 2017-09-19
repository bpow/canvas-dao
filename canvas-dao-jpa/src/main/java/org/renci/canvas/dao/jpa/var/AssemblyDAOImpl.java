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
import org.renci.canvas.dao.var.AssemblyDAO;
import org.renci.canvas.dao.var.model.Assembly;
import org.renci.canvas.dao.var.model.Assembly_;
import org.renci.canvas.dao.var.model.Library;
import org.renci.canvas.dao.var.model.Library_;
import org.renci.canvas.dao.var.model.Sample;
import org.renci.canvas.dao.var.model.Sample_;
import org.renci.canvas.dao.var.model.VariantSet;
import org.renci.canvas.dao.var.model.VariantSet_;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional
@OsgiServiceProvider(classes = { AssemblyDAO.class })
@javax.transaction.Transactional
@Singleton
public class AssemblyDAOImpl extends BaseDAOImpl<Assembly, Integer> implements AssemblyDAO {

    private static final Logger logger = LoggerFactory.getLogger(AssemblyDAOImpl.class);

    public AssemblyDAOImpl() {
        super();
    }

    @Override
    public Class<Assembly> getPersistentClass() {
        return Assembly.class;
    }

    @org.springframework.transaction.annotation.Transactional
    @javax.transaction.Transactional
    @Override
    public void delete(Assembly entity) throws CANVASDAOException {
        logger.debug("ENTERING delete(T)");
        Assembly foundEntity = getEntityManager().find(getPersistentClass(), entity.getId());
        getEntityManager().remove(foundEntity);
    }

    public List<Assembly> findByVariantSetId(Integer variantSetId) throws CANVASDAOException {
        logger.debug("ENTERING findByVariantSetId(Integer)");
        List<Assembly> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<Assembly> crit = critBuilder.createQuery(getPersistentClass());
            Root<Assembly> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Join<Assembly, VariantSet> assemblyVariantSetJoin = root.join(Assembly_.variantSet);
            predicates.add(critBuilder.equal(assemblyVariantSetJoin.get(VariantSet_.id), variantSetId));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<Assembly> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            throw new CANVASDAOException(e);
        }
        return ret;
    }

    @Override
    public List<Assembly> findByLibraryId(Integer libraryId) throws CANVASDAOException {
        logger.debug("ENTERING findByLibraryId(Integer)");
        List<Assembly> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<Assembly> crit = critBuilder.createQuery(getPersistentClass());
            Root<Assembly> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Join<Assembly, Library> assemblyLibraryJoin = root.join(Assembly_.library);
            predicates.add(critBuilder.equal(assemblyLibraryJoin.get(Library_.id), libraryId));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<Assembly> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            throw new CANVASDAOException(e);
        }
        return ret;
    }

    public List<Assembly> findBySampleName(String name) throws CANVASDAOException {
        logger.debug("ENTERING findBySampleName(String)");
        List<Assembly> ret = new ArrayList<>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<Assembly> crit = critBuilder.createQuery(getPersistentClass());
            Root<Assembly> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            Join<Assembly, Library> assemblyLibraryJoin = root.join(Assembly_.library);
            Join<Library, Sample> librarySampleJoin = assemblyLibraryJoin.join(Library_.sample);
            predicates.add(critBuilder.equal(librarySampleJoin.get(Sample_.name), name));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<Assembly> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            throw new CANVASDAOException(e);
        }
        return ret;
    }

}
