package org.tgen.sol.SNP;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.tgen.sol.MappedBaseInfo;
import org.tgen.sol.PositionInfo;

public class SolSNPCaller {

	StrandMode strandMode;
	double callBias;
	Ploidy ploidy;

	public static double[][] TRANSITION_ERROR_PROBABLITY = {{0, 0, 0, 0},
			{0, 0, 0, 0},
			{0, 0, 0, 0},
			{0, 0, 0, 0}};

	final List<MappedBaseInfo> positive_strand_bases = new ArrayList<MappedBaseInfo>();
	final List<MappedBaseInfo> negative_strand_bases = new ArrayList<MappedBaseInfo>();

	public SolSNPCaller(StrandMode sm, double tl, Ploidy pl) {
		strandMode = sm;
		callBias = tl;
		ploidy = pl;
	}

	public SNPCall isSNP(PositionInfo position, char reference) {
		List<MappedBaseInfo> bases;
		switch (strandMode) {
			case None:
				bases = new ArrayList<MappedBaseInfo>();
				for (MappedBaseInfo b : position.mappedBases) {
					bases.add(b);
				}
				return isSNP(bases, reference, callBias);
			case OneStrandAndTotal: {
				bases = new ArrayList<MappedBaseInfo>();

				positive_strand_bases.clear();
				negative_strand_bases.clear();

				for (MappedBaseInfo b : position.mappedBases) {
					if (b.strand)
						positive_strand_bases.add(b);
					else
						negative_strand_bases.add(b);

					bases.add(b);
				}

				SNPCall positive_result = isSNP(positive_strand_bases, reference, callBias);
				SNPCall negative_result = isSNP(negative_strand_bases, reference, callBias);

				SNPCall cons_call = isSNP(bases, reference, callBias);

				cons_call.scores.put("score_ref_positive", positive_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_positive", positive_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_positive", positive_result.scores.get("score_het"));

				cons_call.scores.put("score_ref_negative", negative_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_negative", negative_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_negative", negative_result.scores.get("score_het"));
				cons_call.scores.put("positive_call", (double) positive_result.callType.ordinal());
				cons_call.scores.put("negative_call", (double) negative_result.callType.ordinal());

				int pos_count = 0;
				int neg_count = 0;

				cons_call.misc = "positive_pileup=";
				for (MappedBaseInfo b : positive_strand_bases) {
					pos_count++;
					cons_call.misc += b.nucleotide;
				}

				cons_call.misc += ";negative_pileup=";
				for (MappedBaseInfo b : negative_strand_bases) {
					neg_count++;
					cons_call.misc += b.nucleotide;
				}
				cons_call.misc += ";cov_ratio=";
				if (pos_count == 0 || neg_count == 0)
					cons_call.misc += new Double(0);
				else
					cons_call.misc += new Double(pos_count < neg_count ? pos_count / neg_count : neg_count / pos_count);

				cons_call.misc += ";";

				if (positive_result.callType == cons_call.callType || negative_result.callType == cons_call.callType) {
					return cons_call;
				} else {
					cons_call.callType = CallType.NoCall;
				}
			}
			case NoneWithStrandInfo: {
				bases = new ArrayList<MappedBaseInfo>();
				positive_strand_bases.clear();
				negative_strand_bases.clear();

				for (MappedBaseInfo b : position.mappedBases) {

					if (b.strand)
						positive_strand_bases.add(b);
					else
						negative_strand_bases.add(b);

					bases.add(b);
				}

				SNPCall positive_result = isSNP(positive_strand_bases, reference, callBias);
				SNPCall negative_result = isSNP(negative_strand_bases, reference, callBias);

				SNPCall cons_call = isSNP(bases, reference, callBias);

				cons_call.scores.put("score_ref_positive", positive_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_positive", positive_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_positive", positive_result.scores.get("score_het"));

				cons_call.scores.put("score_ref_negative", negative_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_negative", negative_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_negative", negative_result.scores.get("score_het"));
				cons_call.scores.put("positive_call", (double) positive_result.callType.ordinal());
				cons_call.scores.put("negative_call", (double) negative_result.callType.ordinal());

				int pos_count = 0;
				int neg_count = 0;

				cons_call.misc = "positive_pileup=";
				for (MappedBaseInfo b : positive_strand_bases) {
					pos_count++;
					cons_call.misc += b.nucleotide;
				}

				cons_call.misc += ";negative_pileup=";
				for (MappedBaseInfo b : negative_strand_bases) {
					neg_count++;
					cons_call.misc += b.nucleotide;
				}
				cons_call.misc += ";cov_ratio=";
				if (pos_count == 0 || neg_count == 0)
					cons_call.misc += new Double(0);
				else
					cons_call.misc += new Double(pos_count < neg_count ? pos_count / neg_count : neg_count / pos_count);

				cons_call.misc += ";";

				return cons_call;

			}
			case VariantConsensus: {
				positive_strand_bases.clear();
				negative_strand_bases.clear();

				for (MappedBaseInfo b : position.mappedBases) {
					if (b.strand)
						positive_strand_bases.add(b);
					else
						negative_strand_bases.add(b);
				}

				SNPCall positive_result = isSNP(positive_strand_bases, reference, callBias);
				SNPCall negative_result = isSNP(negative_strand_bases, reference, callBias);

				SNPCall cons_call = new SNPCall();

				cons_call.allele1 = positive_result.allele1;
				cons_call.allele2 = positive_result.allele2;
				cons_call.callType = positive_result.callType;
				cons_call.confidence = (positive_result.confidence + negative_result.confidence) / 2;

				cons_call.scores.put("score_ref_positive", positive_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_positive", positive_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_positive", positive_result.scores.get("score_het"));

				cons_call.scores.put("score_ref_negative", negative_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_negative", negative_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_negative", negative_result.scores.get("score_het"));
				cons_call.scores.put("positive_call", (double) positive_result.callType.ordinal());
				cons_call.scores.put("negative_call", (double) negative_result.callType.ordinal());

				double pos_count = 0;
				double neg_count = 0;

				cons_call.misc = "positive_pileup=";
				for (MappedBaseInfo b : positive_strand_bases) {
					pos_count++;
					cons_call.misc += b.nucleotide;
				}

				cons_call.misc += ";negative_pileup=";
				for (MappedBaseInfo b : negative_strand_bases) {
					neg_count++;
					cons_call.misc += b.nucleotide;
				}
				cons_call.misc += ";cov_ratio=";
				if (pos_count == 0 || neg_count == 0)
					cons_call.misc += 0.0;
				else
					cons_call.misc += new Double(pos_count < neg_count ? pos_count / neg_count : neg_count / pos_count);

				cons_call.misc += ";";

				if (positive_result.callType == CallType.Heterozygote && negative_result.callType == CallType.HomozygoteNonReference
						|| positive_result.callType == CallType.HomozygoteNonReference && negative_result.callType == CallType.Heterozygote) {
					//Indeterminate genotype call, but a SNP call nonetheless
					cons_call.callType = CallType.Heterozygote;

					cons_call.allele1 = reference;

					if (positive_result.allele2 != negative_result.allele2)
						cons_call.allele2 = 'X';
				}
				// Assert that the non-reference alleles in the two strands are equal
				else if (positive_result.callType != negative_result.callType || positive_result.allele2 != negative_result.allele2)
					cons_call.callType = CallType.NoCall;

				if (cons_call.callType == CallType.HomozygoteReference && cons_call.allele2 != cons_call.allele1)
					System.out.println("TWANG");


				return cons_call;
			}
			case GenotypeConsensus: {
				List<MappedBaseInfo> positive_strand_bases = new ArrayList<MappedBaseInfo>();
				List<MappedBaseInfo> negative_strand_bases = new ArrayList<MappedBaseInfo>();

				for (MappedBaseInfo b : position.mappedBases) {

					if (b.strand)
						positive_strand_bases.add(b);
					else
						negative_strand_bases.add(b);
				}

				SNPCall positive_result = isSNP(positive_strand_bases, reference, callBias);
				SNPCall negative_result = isSNP(negative_strand_bases, reference, callBias);

				SNPCall cons_call = new SNPCall();

				cons_call.allele1 = positive_result.allele1;
				cons_call.allele2 = positive_result.allele2;
				cons_call.callType = positive_result.callType;
				cons_call.confidence = (positive_result.confidence + negative_result.confidence) / 2;
				cons_call.scores.put("score_ref_positive", positive_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_positive", positive_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_positive", positive_result.scores.get("score_het"));

				cons_call.scores.put("score_ref_negative", negative_result.scores.get("score_ref"));
				cons_call.scores.put("score_nonref_negative", negative_result.scores.get("score_nonref"));
				cons_call.scores.put("score_het_negative", negative_result.scores.get("score_het"));
				cons_call.scores.put("positive_call", (double) positive_result.callType.ordinal());
				cons_call.scores.put("negative_call", (double) negative_result.callType.ordinal());

				double pos_count = 0;
				double neg_count = 0;

				cons_call.misc = "positive_pileup=";
				for (MappedBaseInfo b : positive_strand_bases) {
					pos_count++;
					cons_call.misc += b.nucleotide;
				}

				cons_call.misc += ";negative_pileup=";
				for (MappedBaseInfo b : negative_strand_bases) {
					neg_count++;
					cons_call.misc += b.nucleotide;
				}
				cons_call.misc += ";cov_ratio=";
				if (pos_count == 0 || neg_count == 0)
					cons_call.misc += new Double(0);
				else
					cons_call.misc += new Double(pos_count < neg_count ? pos_count / neg_count : neg_count / pos_count);

				cons_call.misc += ";";

				// Assert that the non-reference alleles in the two strands are equal
				if (positive_result.callType != negative_result.callType || positive_result.allele2 != negative_result.allele2)
					cons_call.callType = CallType.NoCall;

				return cons_call;
			}
			case NegativeOnly: {
				bases = new ArrayList<MappedBaseInfo>();
				for (MappedBaseInfo b : position.mappedBases) {

					if (!b.strand)
						bases.add(b);
				}
				return isSNP(bases, reference, callBias);
			}
			case PositiveOnly: {
				bases = new ArrayList<MappedBaseInfo>();
				for (MappedBaseInfo b : position.mappedBases) {
					if (b.strand)
						bases.add(b);
				}
				return isSNP(bases, reference, callBias);
			}
			default:
				return null;
		}
	}

	public SNPCall isSNP(List<MappedBaseInfo> mappedBases, char reference, double tolerance) {
		int count = mappedBases.size();
		SNPCall call = new SNPCall();
		call.callType = CallType.NoCall;

		List<MappedBaseInfo> nonreferenceMappedBases = new ArrayList<MappedBaseInfo>();

		char referenceAllele = reference;
		char nonreferenceAllele = ' ';
		char third_nonreference = ' '; // if there is more than one non-reference nucleotide at this locus

		List<Double> refBases = new ArrayList<Double>();
		List<Double> nonrefBases = new ArrayList<Double>();

		for (MappedBaseInfo b : mappedBases) {
			if (b.nucleotide == referenceAllele)
				refBases.add(qualityToProbability(b.quality));
			else {
				if (nonreferenceAllele == ' ')
					nonreferenceAllele = b.nucleotide;
				else if (nonreferenceAllele != b.nucleotide)
					third_nonreference = b.nucleotide; //A third allele was found; Note this so we can resolve which one is the dominant one

				nonrefBases.add(qualityToProbability(b.quality));
				nonreferenceMappedBases.add(b);
			}
		}

		// if there is no evidence of a non-reference allele, we definitely cannot
		// call it a SNP
		if (nonreferenceAllele == ' ') {
			call.allele1 = referenceAllele;
			call.allele2 = referenceAllele;
			call.confidence = 1.0;
			call.callType = CallType.HomozygoteReference;
			return call;
		}

		// if there is evidence of more than one non-reference allele, we need to select one
		if (third_nonreference != ' ') {
			// use the first one we saw as the reference
			// (note: this may not be desirable default behavior
			// if the algorithm is biased towards the reference
			SNPCall nonref_call = isSNP(nonreferenceMappedBases, nonreferenceAllele, 0);

			switch (nonref_call.callType) {
				case HomozygoteNonReference: {
					if (refBases.size() == 0) {
						return nonref_call;
					}
					nonreferenceAllele = third_nonreference;

					break;
				}
				case HomozygoteReference: {
					if (refBases.size() == 0) {
						nonref_call.callType = CallType.HomozygoteNonReference;
						return nonref_call;
					}

					break;
				}
				case Heterozygote: {
					if (refBases.size() == 0) {
						return nonref_call;
					} else {
						nonreferenceAllele = 'X';
					}
				}
				default: {
					if (refBases.size() == 0)
						referenceAllele = 'X';
					nonreferenceAllele = 'X';
				}
			}
		}


		Collections.sort(refBases);
		Collections.sort(nonrefBases);

		// For sample -> reference distribution comparison:
		// Score = The sum of probability of the non-reference bases
		double score_against_reference = sumList(nonrefBases) / count;

		// For sample -> non-reference distribution comparison:
		// Score = The sum of probability of the reference bases
		double score_against_nonreference = sumList(refBases) / count;

		// For sample -> heterozygote distribution comparison:
		// Score = The sum of probability of the 'different' bases, choosing the ones with the best quality
		double score_against_heterozygote = 0;

		switch (ploidy)	 //BUGBUG: Heterozygous score calculated incorrectly if there's only one base mapped
		{
			case Diploid:
				if (refBases.size() > nonrefBases.size()) {
					score_against_heterozygote = sumList(refBases.subList((int) Math.floor(count / 2), refBases.size())) / count;
				} else if (refBases.size() < nonrefBases.size()) {
					score_against_heterozygote = sumList(nonrefBases.subList((int) Math.floor(count / 2), nonrefBases.size())) / count;
				} else
					score_against_heterozygote = Math.abs(sumList(refBases) - sumList(nonrefBases)) / count; // lowest possible score, perfect split (is this right?)
				break;
			case Haploid:
				score_against_heterozygote = 2; //impossible value
				break;
		}

		call.scores.put("score_ref", score_against_reference);
		call.scores.put("score_nonref", score_against_nonreference);
		call.scores.put("score_het", score_against_heterozygote);

		double score_SNP = Math.min(score_against_heterozygote, score_against_nonreference);

		call.scores.put("sens", score_against_reference - score_SNP);

		if (score_SNP <= (score_against_reference + tolerance)) {
			if (score_SNP == score_against_nonreference) {
				call.allele1 = nonreferenceAllele;
				call.allele2 = nonreferenceAllele;
				call.callType = CallType.HomozygoteNonReference;
			} else {
				call.allele1 = referenceAllele;
				call.allele2 = nonreferenceAllele;
				call.callType = CallType.Heterozygote;
				call.scores.put("sens", score_against_reference - score_SNP);
			}

			call.confidence = 1 - score_SNP;

		} else {
			call.allele1 = referenceAllele;
			call.allele2 = referenceAllele;
			call.callType = CallType.HomozygoteReference;
			call.confidence = 1 - score_against_reference;
		}

		return call;
	}

	public double sumList(List<Double> list) {
		double sum = 0;
		for (Double d : list)
			sum += d;
		return sum;
	}

	public double qualityToProbability(int quality) {
		return (1 - Math.pow(10, (-quality / 10)));
	}

}
