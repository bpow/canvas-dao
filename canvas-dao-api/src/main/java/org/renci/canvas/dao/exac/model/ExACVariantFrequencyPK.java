package org.renci.canvas.dao.exac.model;

import java.io.Serializable;

import javax.persistence.Column;
import javax.persistence.Embeddable;

@Embeddable
public class ExACVariantFrequencyPK implements Serializable {

    private static final long serialVersionUID = 4950949109902764215L;

    @Column(name = "loc_var_id")
    private Long locatedVariant;

    @Column(name = "exac_version", length = 10)
    private String version;

    @Column(name = "population", length = 5)
    private String population;

    public ExACVariantFrequencyPK() {
        super();
    }

    public Long getLocatedVariant() {
        return locatedVariant;
    }

    public void setLocatedVariant(Long locatedVariant) {
        this.locatedVariant = locatedVariant;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getPopulation() {
        return population;
    }

    public void setPopulation(String population) {
        this.population = population;
    }

    @Override
    public String toString() {
        return String.format("VariantFrequencyPK [locatedVariant=%s, version=%s, population=%s]", locatedVariant, version, population);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((locatedVariant == null) ? 0 : locatedVariant.hashCode());
        result = prime * result + ((population == null) ? 0 : population.hashCode());
        result = prime * result + ((version == null) ? 0 : version.hashCode());
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
        ExACVariantFrequencyPK other = (ExACVariantFrequencyPK) obj;
        if (locatedVariant == null) {
            if (other.locatedVariant != null)
                return false;
        } else if (!locatedVariant.equals(other.locatedVariant))
            return false;
        if (population == null) {
            if (other.population != null)
                return false;
        } else if (!population.equals(other.population))
            return false;
        if (version == null) {
            if (other.version != null)
                return false;
        } else if (!version.equals(other.version))
            return false;
        return true;
    }

}
