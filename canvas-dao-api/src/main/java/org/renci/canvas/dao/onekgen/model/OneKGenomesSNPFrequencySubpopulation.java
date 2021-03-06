package org.renci.canvas.dao.onekgen.model;

import javax.persistence.Column;
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.MapsId;
import javax.persistence.Table;

import org.renci.canvas.dao.Persistable;
import org.renci.canvas.dao.var.model.LocatedVariant;

@Entity
@Table(schema = "gen1000", name = "snp_freq_subpopulation")
public class OneKGenomesSNPFrequencySubpopulation implements Persistable<OneKGenomesSNPFrequencySubpopulationPK> {

    private static final long serialVersionUID = -3312702331463575620L;

    @EmbeddedId
    private OneKGenomesSNPFrequencySubpopulationPK id;

    @MapsId("locatedVariant")
    @ManyToOne
    @JoinColumn(name = "loc_var_id")
    private LocatedVariant locatedVariant;

    @Column(name = "alt_allele_freq")
    private Double altAlleleFreq;

    @Column(name = "total_allele_count")
    private Integer totalAlleleCount;

    @Column(name = "alt_allele_count")
    private Integer altAlleleCount;

    public OneKGenomesSNPFrequencySubpopulation() {
        super();
    }

    public OneKGenomesSNPFrequencySubpopulationPK getId() {
        return id;
    }

    public void setId(OneKGenomesSNPFrequencySubpopulationPK id) {
        this.id = id;
    }

    public LocatedVariant getLocatedVariant() {
        return locatedVariant;
    }

    public void setLocatedVariant(LocatedVariant locatedVariant) {
        this.locatedVariant = locatedVariant;
    }

    public Double getAltAlleleFreq() {
        return altAlleleFreq;
    }

    public void setAltAlleleFreq(Double altAlleleFreq) {
        this.altAlleleFreq = altAlleleFreq;
    }

    public Integer getTotalAlleleCount() {
        return totalAlleleCount;
    }

    public void setTotalAlleleCount(Integer totalAlleleCount) {
        this.totalAlleleCount = totalAlleleCount;
    }

    public Integer getAltAlleleCount() {
        return altAlleleCount;
    }

    public void setAltAlleleCount(Integer altAlleleCount) {
        this.altAlleleCount = altAlleleCount;
    }

    @Override
    public String toString() {
        return String.format("SNPFrequencySubpopulation [altAlleleFreq=%s, totalAlleleCount=%s, altAlleleCount=%s]", altAlleleFreq,
                totalAlleleCount, altAlleleCount);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((altAlleleCount == null) ? 0 : altAlleleCount.hashCode());
        result = prime * result + ((altAlleleFreq == null) ? 0 : altAlleleFreq.hashCode());
        result = prime * result + ((id == null) ? 0 : id.hashCode());
        result = prime * result + ((totalAlleleCount == null) ? 0 : totalAlleleCount.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        OneKGenomesSNPFrequencySubpopulation other = (OneKGenomesSNPFrequencySubpopulation) obj;
        if (altAlleleCount == null) {
            if (other.altAlleleCount != null)
                return false;
        } else if (!altAlleleCount.equals(other.altAlleleCount))
            return false;
        if (altAlleleFreq == null) {
            if (other.altAlleleFreq != null)
                return false;
        } else if (!altAlleleFreq.equals(other.altAlleleFreq))
            return false;
        if (id == null) {
            if (other.id != null)
                return false;
        } else if (!id.equals(other.id))
            return false;
        if (totalAlleleCount == null) {
            if (other.totalAlleleCount != null)
                return false;
        } else if (!totalAlleleCount.equals(other.totalAlleleCount))
            return false;
        return true;
    }

}
