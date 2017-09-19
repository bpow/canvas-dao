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
import org.renci.canvas.dao.var.LibraryDAO;
import org.renci.canvas.dao.var.model.Library;
import org.renci.canvas.dao.var.model.Library_;
import org.renci.canvas.dao.var.model.Sample;
import org.renci.canvas.dao.var.model.Sample_;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import org.springframework.transaction.annotation.Transactional;

@Component
@org.springframework.transaction.annotation.Transactional
@OsgiServiceProvider(classes = { LibraryDAO.class })
@javax.transaction.Transactional
@Singleton
public class LibraryDAOImpl extends BaseDAOImpl<Library, Integer> implements LibraryDAO {

    private static final Logger logger = LoggerFactory.getLogger(LibraryDAOImpl.class);

    public LibraryDAOImpl() {
        super();
    }

    @Override
    public Class<Library> getPersistentClass() {
        return Library.class;
    }

    @org.springframework.transaction.annotation.Transactional
    @javax.transaction.Transactional
    @Override
    public void delete(Library entity) throws CANVASDAOException {
        logger.debug("ENTERING delete(T)");
        Library foundEntity = getEntityManager().find(getPersistentClass(), entity.getId());
        getEntityManager().remove(foundEntity);
    }

    @Override
    @Transactional(readOnly = true)
    public List<Library> findByNameAndSampleId(String name, Integer sampleId) throws CANVASDAOException {
        logger.debug("ENTERING findByNameAndSampleId(String, Integer)");
        List<Library> ret = new ArrayList<Library>();
        try {
            CriteriaBuilder critBuilder = getEntityManager().getCriteriaBuilder();
            CriteriaQuery<Library> crit = critBuilder.createQuery(getPersistentClass());
            Root<Library> root = crit.from(getPersistentClass());
            List<Predicate> predicates = new ArrayList<Predicate>();
            predicates.add(critBuilder.equal(root.get(Library_.name), name));
            Join<Library, Sample> librarySampleJoin = root.join(Library_.sample);
            predicates.add(critBuilder.equal(librarySampleJoin.get(Sample_.id), sampleId));
            crit.where(predicates.toArray(new Predicate[predicates.size()]));
            TypedQuery<Library> query = getEntityManager().createQuery(crit);
            ret.addAll(query.getResultList());
        } catch (Exception e) {
            throw new CANVASDAOException(e);
        }
        return ret;
    }

}
