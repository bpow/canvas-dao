package org.renci.canvas.dao.onekgen.model;

import java.io.Serializable;

import javax.persistence.Column;
import javax.persistence.Embeddable;

@Embeddable
public class OneKGenomesSNPPopulationMaxFrequencyPK implements Serializable {

    private static final long serialVersionUID = -3681815447860895929L;

    @Column(name = "loc_var_id")
    private Long locatedVariant;

    @Column(name = "gen1000_version")
    private Integer version;

    public OneKGenomesSNPPopulationMaxFrequencyPK() {
        super();
    }

    public OneKGenomesSNPPopulationMaxFrequencyPK(Long locatedVariant, Integer version) {
        super();
        this.locatedVariant = locatedVariant;
        this.version = version;
    }

    public Long getLocatedVariant() {
        return locatedVariant;
    }

    public void setLocatedVariant(Long locatedVariant) {
        this.locatedVariant = locatedVariant;
    }

    public Integer getVersion() {
        return version;
    }

    public void setVersion(Integer version) {
        this.version = version;
    }

    @Override
    public String toString() {
        return String.format("SNPFrequencySubpopulationPK [locatedVariant=%s, gen1000Version=%s]", locatedVariant, version);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((version == null) ? 0 : version.hashCode());
        result = prime * result + ((locatedVariant == null) ? 0 : locatedVariant.hashCode());
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
        OneKGenomesSNPPopulationMaxFrequencyPK other = (OneKGenomesSNPPopulationMaxFrequencyPK) obj;
        if (version == null) {
            if (other.version != null)
                return false;
        } else if (!version.equals(other.version))
            return false;
        if (locatedVariant == null) {
            if (other.locatedVariant != null)
                return false;
        } else if (!locatedVariant.equals(other.locatedVariant))
            return false;
        return true;
    }

}
