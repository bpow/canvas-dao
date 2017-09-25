package org.renci.canvas.dao.commons;

import java.util.ArrayList;
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

                locatedVariant = createSNP(genomeRef, genomeRefSeq, snp, ref, alt, oneBasedStart);

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
                locatedVariant.setRef("");
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

    public static LocatedVariant createSNP(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq, VariantType snpVariantType, String ref,
            String alt, Integer position) {
        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);
        locatedVariant.setSeq(alt);
        locatedVariant.setRef(ref);
        locatedVariant.setPosition(position);
        locatedVariant.setVariantType(snpVariantType);
        locatedVariant.setEndPosition(position + locatedVariant.getRef().length());
        return locatedVariant;
    }

    public static LocatedVariant createDeletion(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq, VariantType delVariantType, String ref,
            String alt, Integer position) {

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

        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);

        locatedVariant.setVariantType(delVariantType);
        locatedVariant.setPosition(position);

        List<Character> referenceChars = ref.chars().mapToObj(c -> (char) c).collect(Collectors.toList());
        List<Character> alternateChars = alt.chars().mapToObj(c -> (char) c).collect(Collectors.toList());

        if (referenceChars.size() > 1 && alternateChars.size() > 1) {

            StringBuilder frontChars2Remove = new StringBuilder();
            StringBuilder backChars2Remove = new StringBuilder();

            for (int i = 0; i < referenceChars.size(); ++i) {
                if (i == alternateChars.size() || referenceChars.get(i) != alternateChars.get(i)) {
                    break;
                }
                frontChars2Remove.append(referenceChars.get(i));
            }

            if (CollectionUtils.isNotEmpty(referenceChars) && CollectionUtils.isNotEmpty(alternateChars)) {

                String tmpRef = referenceChars.stream().map(a -> a.toString()).collect(Collectors.joining());
                tmpRef = tmpRef.replaceFirst(frontChars2Remove.toString(), "");
                referenceChars = tmpRef.chars().mapToObj(c -> (char) c).collect(Collectors.toList());

                String tmpAlt = alternateChars.stream().map(a -> a.toString()).collect(Collectors.joining());
                tmpAlt = tmpAlt.replaceFirst(frontChars2Remove.toString(), "");
                alternateChars = tmpAlt.chars().mapToObj(c -> (char) c).collect(Collectors.toList());

                Collections.reverse(referenceChars);
                Collections.reverse(alternateChars);
                for (int i = 0; i < referenceChars.size(); ++i) {
                    if (i == alternateChars.size() || referenceChars.get(i) != alternateChars.get(i)) {
                        break;
                    }
                    backChars2Remove.append(referenceChars.get(i));
                }
            }

            if (frontChars2Remove.length() > 0) {
                ref = ref.replaceFirst(frontChars2Remove.toString(), "");
                alt = alt.replaceFirst(frontChars2Remove.toString(), "");
            }

            if (backChars2Remove.length() > 0) {
                backChars2Remove.reverse();
                ref = StringUtils.removeEnd(ref, backChars2Remove.toString());
                alt = StringUtils.removeEnd(alt, backChars2Remove.toString());
            }

            locatedVariant.setPosition(position + (frontChars2Remove.length() > 0 ? frontChars2Remove.length() : 0));

        }

        locatedVariant.setEndPosition(locatedVariant.getPosition() + ref.length());
        locatedVariant.setRef(ref);
        locatedVariant.setSeq(ref);

        return locatedVariant;
    }

    public static LocatedVariant createInsertion(GenomeRef genomeRef, GenomeRefSeq genomeRefSeq, VariantType insVariantType, String ref,
            String alt, Integer position) {

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

        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);

        locatedVariant.setVariantType(insVariantType);
        locatedVariant.setPosition(position - 1);

        List<Character> referenceChars = ref.chars().mapToObj(c -> (char) c).collect(Collectors.toList());
        List<Character> alternateChars = alt.chars().mapToObj(c -> (char) c).collect(Collectors.toList());

        if (referenceChars.size() > 1 && alternateChars.size() > 1) {

            StringBuilder frontChars2Remove = new StringBuilder();
            StringBuilder backChars2Remove = new StringBuilder();

            List<Character> charsToRemoveFromRef = new ArrayList<>();

            for (int i = 0; i < referenceChars.size(); ++i) {
                if (referenceChars.get(i) != alternateChars.get(i)) {
                    break;
                }
                charsToRemoveFromRef.add(referenceChars.get(i));
                frontChars2Remove.append(referenceChars.get(i));
            }

            for (int i = 0; i < charsToRemoveFromRef.size(); ++i) {
                referenceChars.remove(charsToRemoveFromRef.get(i));
            }

            if (CollectionUtils.isNotEmpty(referenceChars)) {
                Collections.reverse(referenceChars);
                Collections.reverse(alternateChars);
                for (int i = 0; i < referenceChars.size(); ++i) {
                    if (referenceChars.get(i) != alternateChars.get(i)) {
                        break;
                    }
                    backChars2Remove.append(referenceChars.get(i));
                }
            }

            if (frontChars2Remove.length() > 0) {
                ref = ref.replaceFirst(frontChars2Remove.toString(), "");
                alt = alt.replaceFirst(frontChars2Remove.toString(), "");
            }

            if (backChars2Remove.length() > 0) {
                backChars2Remove.reverse();
                ref = StringUtils.removeEnd(ref, backChars2Remove.toString());
                alt = StringUtils.removeEnd(alt, backChars2Remove.toString());
            }

            locatedVariant.setPosition(position - 1 + (frontChars2Remove.length() > 0 ? frontChars2Remove.length() : 0));
        }

        locatedVariant.setEndPosition(locatedVariant.getPosition() + 1);
        locatedVariant.setRef("");
        locatedVariant.setSeq(alt);
        return locatedVariant;

    }

}
