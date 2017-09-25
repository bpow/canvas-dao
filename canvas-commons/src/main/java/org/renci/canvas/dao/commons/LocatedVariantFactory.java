package org.renci.canvas.dao.commons;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.renci.canvas.dao.ref.model.GenomeRef;
import org.renci.canvas.dao.ref.model.GenomeRefSeq;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.renci.canvas.dao.var.model.VariantType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class LocatedVariantFactory {

    private static final Logger logger = LoggerFactory.getLogger(LocatedVariantFactory.class);
    private static final List<VariantType> dummyVariantTypes = Arrays.asList(
            new VariantType("snp"),
            new VariantType("del"),
            new VariantType("ins"),
            new VariantType("sub"),
            new VariantType("ref"));

    private static int sharedPrefixLength(String a, String b) {
        int maxShared = Math.min(a.length(), b.length());
        for (int i = 0; i < maxShared; i++) {
            if (a.charAt(i) != b.charAt(i)) return i;
        }
        return maxShared;
    }

    public static LocatedVariant create(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq,
                                        VariantContext variantContext, Allele altAllele,
                                        List<VariantType> allVariantTypes) {
        return create(genomeRef, genomeRefSeq, variantContext.getReference().getDisplayString(),
                altAllele.getDisplayString(), variantContext.getStart(), allVariantTypes);
    }

    public static LocatedVariant create(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq,
                                        String ref, String alt, int oneBasedStart,
                                        List<VariantType> allVariantTypes) {

        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);

        VariantType snp = allVariantTypes.stream().filter(a -> a.getId().equals("snp")).findAny().get();
        VariantType del = allVariantTypes.stream().filter(a -> a.getId().equals("del")).findAny().get();
        VariantType ins = allVariantTypes.stream().filter(a -> a.getId().equals("ins")).findAny().get();
        VariantType sub = allVariantTypes.stream().filter(a -> a.getId().equals("sub")).findAny().get();

        try {

            if (StringUtils.isEmpty(ref)) {
                logger.error("ref is empty");
                return null;
            }

            if (StringUtils.isEmpty(alt)) {
                logger.error("alt is empty");
                return null;
            }

            if (ref.equals(alt)) {
                logger.error("ref and alt are equal");
                return null;
            }

            int prefixLength = sharedPrefixLength(ref, alt);

            String reversedRefWithoutPrefix = new StringBuffer(ref.substring(prefixLength)).reverse().toString();
            String reversedAltWithoutPrefix = new StringBuffer(ref.substring(prefixLength)).reverse().toString();

            int suffixLength = sharedPrefixLength(reversedRefWithoutPrefix, reversedAltWithoutPrefix);

            ref = ref.substring(prefixLength, ref.length()-suffixLength);
            alt = alt.substring(prefixLength, alt.length()-suffixLength);

            oneBasedStart = oneBasedStart + prefixLength;

            if (ref.length() == 1 && alt.length() == 1) { // SNP (or ref...)
                locatedVariant.setVariantType(snp);
                locatedVariant.setPosition(oneBasedStart);
                locatedVariant.setEndPosition(locatedVariant.getPosition() + ref.length());
                locatedVariant.setRef(ref);
                locatedVariant.setSeq(alt);
            } else if (alt.length() == 0) { // deletion
                locatedVariant.setVariantType(del);
                locatedVariant.setPosition(oneBasedStart);
                locatedVariant.setEndPosition(locatedVariant.getPosition() + ref.length());
                locatedVariant.setRef(ref);
                locatedVariant.setSeq(ref);
            } else if (ref.length() == 0) { // insertion
                locatedVariant.setVariantType(ins);
                locatedVariant.setPosition(oneBasedStart-1);
                locatedVariant.setEndPosition(locatedVariant.getPosition() + 1);
                locatedVariant.setRef(ref);
                locatedVariant.setSeq(alt);
            } else if (ref.length() == alt.length()) { // sub (delins)
                locatedVariant.setVariantType(sub);
                locatedVariant.setPosition((oneBasedStart-1));
                locatedVariant.setRef(ref);
                locatedVariant.setSeq(alt);
                locatedVariant.setEndPosition(locatedVariant.getPosition() + ref.length());
            }

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

        return locatedVariant;
    }

    /**
     * @deprecated Use {@link #create(GenomeRef, GenomeRefSeq, String, String, int, List)} instead
     */
    @Deprecated public static LocatedVariant createSNP(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq, VariantType snpVariantType, String ref,
            String alt, Integer position) {
        LocatedVariant locatedVariant = create(genomeRef, genomeRefSeq, ref, alt, position, dummyVariantTypes);
        // check that we got the expected variant type
        if (!snpVariantType.getId().equals(locatedVariant.getVariantType().getId())) {
            logger.error(String.format("Expected to get a variant with type '%s', got variant with type '%s' instead",
                    snpVariantType.getId(), locatedVariant.getVariantType().getId()));
            return null;
        }
        // use the passed-in VariantType so persistence will work
        locatedVariant.setVariantType(snpVariantType);
        return locatedVariant;
    }

    /**
     * @deprecated Use {@link #create(GenomeRef, GenomeRefSeq, String, String, int, List)} instead
     */
    @Deprecated public static LocatedVariant createDeletion(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq, VariantType delVariantType, String ref,
            String alt, Integer position) {
        LocatedVariant locatedVariant = create(genomeRef, genomeRefSeq, ref, alt, position, dummyVariantTypes);
        // check that we got the expected variant type
        if (!delVariantType.getId().equals(locatedVariant.getVariantType().getId())) {
            logger.error(String.format("Expected to get a variant with type '%s', got variant with type '%s' instead",
                    delVariantType.getId(), locatedVariant.getVariantType().getId()));
            return null;
        }
        // use the passed-in VariantType so persistence will work
        locatedVariant.setVariantType(delVariantType);
        return locatedVariant;
    }

    /**
     * @deprecated Use {@link #create(GenomeRef, GenomeRefSeq, String, String, int, List)} instead
     */
    @Deprecated public static LocatedVariant createInsertion(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq, VariantType insVariantType, String ref,
            String alt, Integer position) {
        LocatedVariant locatedVariant = create(genomeRef, genomeRefSeq, ref, alt, position, dummyVariantTypes);
        // check that we got the expected variant type
        if (!insVariantType.getId().equals(locatedVariant.getVariantType().getId())) {
            logger.error(String.format("Expected to get a variant with type '%s', got variant with type '%s' instead",
                    insVariantType.getId(), locatedVariant.getVariantType().getId()));
            return null;
        }
        // use the passed-in VariantType so persistence will work
        locatedVariant.setVariantType(insVariantType);
        return locatedVariant;
    }

}
