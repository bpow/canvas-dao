package org.renci.canvas.dao.jpa.refseq;

import java.util.List;

import javax.inject.Singleton;
import javax.persistence.TypedQuery;

import org.ops4j.pax.cdi.api.OsgiServiceProvider;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.renci.canvas.dao.refseq.Variants_48_2_DAO;
import org.renci.canvas.dao.refseq.model.Variants_48_2;
import org.renci.canvas.dao.refseq.model.Variants_48_2PK;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional(readOnly = true)
@OsgiServiceProvider(classes = { Variants_48_2_DAO.class })
@javax.transaction.Transactional(javax.transaction.Transactional.TxType.SUPPORTS)
@Singleton
public class Variants_48_2_DAOImpl extends BaseDAOImpl<Variants_48_2, Variants_48_2PK> implements Variants_48_2_DAO {

    private static final Logger logger = LoggerFactory.getLogger(Variants_48_2_DAOImpl.class);

    public Variants_48_2_DAOImpl() {
        super();
    }

    @Override
    public Class<Variants_48_2> getPersistentClass() {
        return Variants_48_2.class;
    }

    @Override
    public List<Variants_48_2> findByLocatedVariantId(Long id) throws CANVASDAOException {
        logger.debug("ENTERING findByName()");
        TypedQuery<Variants_48_2> query = getEntityManager().createNamedQuery("Variants_48_2.findByLocatedVariantId", Variants_48_2.class);
        query.setParameter("LocatedVariantId", id);
        List<Variants_48_2> ret = query.getResultList();
        return ret;
    }

}
