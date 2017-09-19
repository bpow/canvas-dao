package org.renci.canvas.dao.jpa.dbsnp;

import javax.inject.Singleton;

import org.ops4j.pax.cdi.api.OsgiServiceProvider;
import org.renci.canvas.dao.dbsnp.SNPAlleleDAO;
import org.renci.canvas.dao.dbsnp.model.SNPAllele;
import org.renci.canvas.dao.dbsnp.model.SNPAllelePK;
import org.renci.canvas.dao.jpa.BaseDAOImpl;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
@org.springframework.transaction.annotation.Transactional(readOnly = true)
@OsgiServiceProvider(classes = { SNPAlleleDAO.class })
@javax.transaction.Transactional(javax.transaction.Transactional.TxType.SUPPORTS)
@Singleton
public class SNPAlleleDAOImpl extends BaseDAOImpl<SNPAllele, SNPAllelePK> implements SNPAlleleDAO {

    private static final Logger logger = LoggerFactory.getLogger(SNPAlleleDAOImpl.class);

    public SNPAlleleDAOImpl() {
        super();
    }

    @Override
    public Class<SNPAllele> getPersistentClass() {
        return SNPAllele.class;
    }

}
