/**
 *
 */
package org.tgen.sol.SNP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

import org.tgen.sol.*;


enum StrandMode {
	GenotypeConsensus,
	VariantConsensus,
	PositiveOnly,
	NegativeOnly,
	None,
	NoneWithStrandInfo,
	OneStrandAndTotal
}

enum OutputFormat {
	GFF,
	VCF
}

enum Ploidy {
	Haploid,
	Diploid
}

/**
 * @author achristoforides
 */
public class SolSNP extends CommandLineProgram {

	@Usage
	public final String USAGE = getStandardUsagePreamble() + "SolSNP: Uses a modified Kolmogorov–Smirnov test to produce SNP variant calls from a SAM/BAM alignment file and a reference.";

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM/BAM file. Needs to be sorted by coordinate.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file.")
	public File OUTPUT;

	@Option(doc = "Known calls file.", optional = true)
	public File KNOWN_CALLS;

	@Option(doc = "Strand Mode", optional = true)
	public StrandMode STRAND_MODE = StrandMode.VariantConsensus;

	@Option(doc = "Set to true to create a summary metrics directory upon completion of analysis.", optional = true)
	public Boolean SUMMARY = false;

	@Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME,
			doc = "Reference Sequence File")
	public File REFERENCE_SEQUENCE;

	@Option(doc = "Minimum confidence score allowed for calls. ")
	public double FILTER;

	@Option(doc = "Additional score bias towards a variant call. (Range: 0.0 - 1.0)", optional = true)
	public double CALL_BIAS = 0;

	@Option(doc = "Minimum base quality", optional = true)
	public int MINIMUM_BASE_QUALITY = 1; // By default, only bases with quality zero are trimmed

	@Option(doc = "Minimum mapping quality", optional = true)
	public int MINIMUM_MAPQ = 1; // By default, only reads with mapping quality zero are trimmed

	@Option(doc = "Output Format", optional = true)
	public OutputFormat OUTPUT_FORMAT = OutputFormat.GFF;

	@Option(doc = "Ploidy", optional = true)
	public Ploidy PLOIDY = Ploidy.Diploid;

	@Option(doc = "Minimum coverage", optional = true)
	public static int MINIMUM_COVERAGE = 3;

	@Option(doc = "Region (syntax: sequence_name,start,end)", optional = true)
	public static String REGION = "";

	@Option(doc = "Set to true to include non-variant genotype calls to the output file.", optional = true)
	public Boolean OUTPUT_REFERENCE_MATCHES = false;

	@Option(doc = "Maximum genomic distance between paired reads.", optional = true)
	public static int MAX_MATE_DISTANCE = Integer.MAX_VALUE;

	@Option(doc = "Minimum genomic distance between paired reads.", optional = true)
	public static int MIN_MATE_DISTANCE = 0;

	/* (non-Javadoc)
		 * @see net.sf.picard.cmdline.CommandLineProgram#doWork()
		 */

	@Override
	protected int doWork() {
		IoUtil.assertFileIsReadable(INPUT);
		ReferenceSequenceFile reference_file = null;
		if (REFERENCE_SEQUENCE != null) {
			IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
			reference_file = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

		}
		PrintWriter out;
		PrintWriter false_negatives_file = null;

		FileReader known_calls_filereader;
		BufferedReader known_calls_bufferedreader;
		HashMap<String, HashMap<Integer, String>> known_calls = new HashMap<String, HashMap<Integer, String>>();


		if (KNOWN_CALLS != null) {
			try {
				known_calls_filereader = new FileReader(KNOWN_CALLS);
				known_calls_bufferedreader = new BufferedReader(known_calls_filereader);
				false_negatives_file = new PrintWriter(OUTPUT + ".falsenegatives");

				String currentRecord;

				while ((currentRecord = known_calls_bufferedreader.readLine()) != null) {
					String[] fields = (currentRecord.split("\t"));
					if (fields.length == 4) {
						if (!known_calls.containsKey(fields[1]))
							known_calls.put(fields[1], new HashMap<Integer, String>());

						known_calls.get(fields[1]).put(Integer.parseInt(fields[2]), fields[3]);
					} else
						throw new PicardException("Invalid 'known calls' format.");
				}
			}

			catch (Exception e) {
				throw new PicardException(e.getMessage());
			}

		}

		if (OUTPUT != null) {
			IoUtil.assertFileIsWritable(OUTPUT);
			try {
				out = new PrintWriter(OUTPUT);

			}
			catch (FileNotFoundException e) {
				// we already asserted this so we should not get here
				throw new PicardException("Unexpected exception", e);
			}
		} else {
			out = new PrintWriter(System.out);
		}

		PrintWriter summary_file = null;
		PrintWriter strand_coverage_file = null;
		PrintWriter validation_file = null;

		PrintWriter validation_coverage_file = null;

		if (SUMMARY) {
			String summarydir_name = OUTPUT + "_summary" + File.separator;
			new File(summarydir_name).mkdir();
			try {
				summary_file = new PrintWriter(summarydir_name + "summary.txt");
				strand_coverage_file = new PrintWriter(summarydir_name + "strandcoveragetable.txt");
				validation_file = new PrintWriter(summarydir_name + "validation.txt");
				validation_coverage_file = new PrintWriter(summarydir_name + "validation_coverage.txt");
			}
			catch (FileNotFoundException e) {
				// we already asserted this so we should not get here
				throw new PicardException("Could not create metrics files", e);
			}

		}

		SamPositionIterator position_iterator;

		if (REGION.equals("")) {
			position_iterator = new SamPositionIterator(INPUT, MINIMUM_COVERAGE, MINIMUM_MAPQ, MIN_MATE_DISTANCE, MAX_MATE_DISTANCE, MINIMUM_BASE_QUALITY);
		} else {
			//split the region string into sequence name/beginning/end
			String[] tokens = REGION.split(",");
			int start = 0, end = 0;

			try {
				start = Integer.parseInt(tokens[1]);
				end = Integer.parseInt(tokens[2]);
			}
			catch (Exception e) {
				System.err.println("Error parsing definition string");
			}
			position_iterator = new SamPositionIterator(INPUT, MINIMUM_COVERAGE, MINIMUM_MAPQ, MIN_MATE_DISTANCE, MAX_MATE_DISTANCE, MINIMUM_BASE_QUALITY, tokens[0], start, end);

		}
		//SAMSequenceDictionary reference_dictionary = reference_file.getSequenceDictionary();


		int snp_count = 0;
		int transitions = 0;
		int transversions = 0;

		double balance_ref = 0.0;
		double balance_het = 0.0;
		double balance_hom = 0.0;
		int n_ref = 0;
		int n_het = 0;
		int n_hom = 0;

		int mismatch_total = 0;
		int max_coverage = 0;

		SolSNPCaller s = new SolSNPCaller(STRAND_MODE, CALL_BIAS, PLOIDY);

		int coverage;
		int previous_position = 0;
		int calls_at_coverage[] = new int[65536];
		int wrongcalls_at_coverage[] = new int[65536];
		int missedcalls_at_coverage[] = new int[65536];

		int strand_coverage[][] = new int[50][50];

		final HashMap<String, Integer> mismatch_counts = new HashMap<String, Integer>();
		final Map<SNPCallPair, Long> calltable = new HashMap<SNPCallPair, Long>();
		final Map<SNPCallPair, Long> mismatchtable = new HashMap<SNPCallPair, Long>();
		final Map<SNPCallPair, TreeMap<Integer, Long>> calltable_validation = new HashMap<SNPCallPair, TreeMap<Integer, Long>>();
		final SNPCallPair lookup = new SNPCallPair(CallType.Uncallable, CallType.Uncallable);

		Set<Character> nuc_values = new TreeSet<Character>();

		while (position_iterator.nextSequence()) {
			ReferenceSequence reference = reference_file.nextSequence();

			if (reference == null) {
				System.out.println("Warning: Input file contains records mapped on a sequence that is not in the reference file: " + position_iterator.getCurrentSequence());
				break;
			}
			String ref_name = reference.getName();

			if (!position_iterator.getCurrentSequence().equals(ref_name))
				continue;

			System.out.println("Processing reads on sequence " + ref_name);

			byte[] ref_array = reference.getBases();

			reference = null;

			HashMap<Integer, String> ref_known_calls = known_calls.get(ref_name);

			if (ref_known_calls == null)
				ref_known_calls = new HashMap<Integer, String>();

			// Iterate through loci with enough coverage
			while (position_iterator.hasMoreElements()) {
				PositionInfo p = position_iterator.nextElement();

				int curp;

				if (p != null)
					curp = p.position;
				else
					curp = ref_array.length;

				if (p != null) {
					char reference_nucleotide = (char) ref_array[p.position - 1];
					Iterator<MappedBaseInfo> i = p.mappedBases.iterator();
					while (i.hasNext()) {
						MappedBaseInfo b = i.next();

						if (b.nucleotide != reference_nucleotide) {
							mismatch_total++;
							String mismatch_lookup = "";
							mismatch_lookup += reference_nucleotide;
							mismatch_lookup += b.nucleotide; //new String(new char[] {reference_nucleotide, b.nucleotide , });

							Integer n = null;

							n = mismatch_counts.get(mismatch_lookup);

							if (n == null) {
								n = new Integer(0);
							}
							n = n + 1;
							mismatch_counts.put(mismatch_lookup, n);
						}
					}

					if (p.mappedBases.size() < MINIMUM_COVERAGE)
						continue;
				}


				for (int nc = previous_position + 1; nc < curp; nc++) {
					//System.out.println(nc);
					Character reference_nucleotide = (char) ref_array[curp - 1];
					//uncallable_count += (p.position - previous_position - 1);
					String known_call = ref_known_calls.get(nc);
					CallType known_call_type = CallType.Unknown;

					if (known_call != null) // call
					{
						if (known_call.charAt(0) == known_call.charAt(1)) {
							if (known_call.charAt(0) == reference_nucleotide)
								known_call_type = CallType.HomozygoteReference;

							else
								known_call_type = CallType.HomozygoteNonReference;

						} else
							known_call_type = CallType.Heterozygote;

					}


					lookup.call1 = CallType.Uncallable;
					lookup.call2 = known_call_type;

					Long count = calltable.get(lookup);
					if (count == null) {
						count = new Long(1);
						calltable.put((SNPCallPair) lookup.clone(), count);
					} else {
						count = count + 1;
						calltable.put(lookup, count);
					}


				}
				previous_position = curp;

				if (p == null)
					break;

				Character reference_nucleotide = (char) ref_array[p.position - 1];

				String known_call = ref_known_calls.get(p.position);
				coverage = p.mappedBases.size();
				if (coverage > max_coverage)
					max_coverage = coverage;

				CallType known_call_type = CallType.Unknown;

				if (known_call != null) {
					if (known_call.charAt(0) == known_call.charAt(1)) {
						if (known_call.charAt(0) == reference_nucleotide) {
							known_call_type = CallType.HomozygoteReference;
							balance_ref += CalculateAlleleBalance(p.mappedBases, reference_nucleotide);
							n_ref++;
						} else {
							known_call_type = CallType.HomozygoteNonReference;
							balance_hom += CalculateAlleleBalance(p.mappedBases, reference_nucleotide);
							n_hom++;
						}
					} else {
						known_call_type = CallType.Heterozygote;
						balance_het += CalculateAlleleBalance(p.mappedBases, reference_nucleotide);
						n_het++;
					}
				}

				//Attempt to call SNP at this locus
				SNPCall SNP = s.isSNP(p, reference_nucleotide);

				int pos_count = 0;
				int neg_count = 0;
				for (MappedBaseInfo m : p.mappedBases) {
					if (m.strand)
						pos_count++;
					else
						neg_count++;
				}
				if (pos_count < 50 && neg_count < 50)
					strand_coverage[pos_count][neg_count]++;

				//Apply low-confidence filter
				if (SNP.confidence < FILTER)
					SNP.callType = CallType.NoCall;

				lookup.call1 = SNP.callType;
				lookup.call2 = known_call_type;


				Long count = calltable.get(lookup);
				if (count == null) {
					count = new Long(1);
					calltable.put((SNPCallPair) lookup.clone(), count);
				} else {
					count = count + 1;
					calltable.put(lookup, count);
				}

				TreeMap<Integer, Long> map = calltable_validation.get(lookup);
				if (map == null) {
					map = new TreeMap<Integer, Long>();
					calltable_validation.put(lookup, map);
				}
				count = map.get(coverage);
				if (count == null) {
					count = new Long(1);
					map.put(coverage, count);
				} else {
					count = count + 1;
					map.put(coverage, count);
				}

				if ((known_call_type != CallType.Unknown && (SNP.callType != known_call_type || !(known_call.equals(SNP.toString()) || known_call.equals(SNP.toString2()))))) {
					Long mismatchcount = mismatchtable.get(lookup);
					if (mismatchcount == null) {
						mismatchcount = new Long(1);
						mismatchtable.put(lookup, mismatchcount);
					} else {
						mismatchcount = mismatchcount + 1;
						mismatchtable.put(lookup, mismatchcount);
					}

					wrongcalls_at_coverage[coverage]++;
				}

				//count validation calls at coverage
				if (known_call_type != CallType.Unknown && known_call_type != CallType.NoCall)
					calls_at_coverage[coverage]++;

				if (SNP.callType == CallType.HomozygoteNonReference || SNP.callType == CallType.Heterozygote || OUTPUT_REFERENCE_MATCHES && (SNP.callType == CallType.HomozygoteReference || SNP.callType == CallType.NoCall)) {
					// transition/transversion count

					nuc_values.clear();

					nuc_values.add(reference_nucleotide);
					nuc_values.add(SNP.allele1);
					nuc_values.add(SNP.allele2);

					if (nuc_values.contains('G') && nuc_values.contains('A'))
						transitions++;
					else if (nuc_values.contains('C') && nuc_values.contains('T'))
						transitions++;
					else
						transversions++;

					switch (OUTPUT_FORMAT) {
						case GFF: {
							out.printf("snp_%s_%d\tks-snp-call\tsnp\t%d\t%d\t%5.6f\t.\t.\tcall=%s;i=%s;ref=%s;",
									p.sequenceName, ++snp_count,
									p.position, p.position, SNP.confidence,
									SNP.toString(), p.sequenceName,
									(char) ref_array[p.position - 1]);

							out.print(SNP.misc);

							for (String score : SNP.scores.keySet()) {
								out.printf("%s=%5.9f;", score, SNP.scores.get(score));
							}

							out.print("pileup=");

							for (MappedBaseInfo b : p.mappedBases)
								out.print(b.nucleotide);

							out.println();
						}
						break;
						case VCF: {
							double phred_score = (SNP.confidence == 1 ? 30 : -10 * Math.log10(1 - SNP.confidence));
							out.printf("%s\t%d\t.\t%s\t%s\t%5.1f\t0\tCL=%s;",
									p.sequenceName,
									p.position,
									(char) ref_array[p.position - 1],
									(SNP.allele2 == (char) ref_array[p.position - 1] ? '.' : SNP.allele2),
									phred_score,
									SNP.toString()
							);

							if (known_call != null)
								out.printf("KC=%s;", known_call);

							out.print("PL=");

							for (MappedBaseInfo b : p.mappedBases)
								if (b.strand)
									out.print(b.nucleotide);
								else
									out.print(Character.toLowerCase(b.nucleotide));

							out.println();
						}
					}
				} else if (known_call_type != CallType.Unknown && known_call_type != CallType.HomozygoteReference && (SNP.callType == CallType.NoCall || SNP.callType == CallType.HomozygoteReference)) {
					missedcalls_at_coverage[coverage]++;

					false_negatives_file.printf("snp_%s_%d\tks-snp-call\tsnp\t%d\t%d\t%5.6f\t.\t.\tcall=%s;known_call=%s;i=%s;ref=%s;",
							p.sequenceName, 0,
							p.position, p.position, SNP.confidence,
							SNP.toString(), known_call, p.sequenceName, (char) ref_array[p.position - 1]);

					false_negatives_file.print(SNP.misc);

					for (String score : SNP.scores.keySet()) {
						false_negatives_file.printf("%s=%5.9f;", score, SNP.scores.get(score));
					}

					false_negatives_file.print("pileup=");

					for (MappedBaseInfo b : p.mappedBases)
						false_negatives_file.print(b.nucleotide);

					boolean firstqual = true;

					false_negatives_file.print(";quals=");

					for (MappedBaseInfo b : p.mappedBases) {
						if (!firstqual)
							false_negatives_file.print(",");

						false_negatives_file.print(b.quality);
						firstqual = false;
					}

					false_negatives_file.println();
				}

			}
		}
		out.close();

		if (KNOWN_CALLS != null)
			false_negatives_file.close();

		// write summary metrics
		if (SUMMARY) {
			//summary

			//Ti/Tv
			summary_file.write("Transition count\t" + transitions + "\n");
			summary_file.write("Transversion count\t" + transversions + "\n");
			;

			//allele balance
			summary_file.write("\nAllelic Balance on known calls\n--\n");
			summary_file.write("Ref\t " + balance_ref / n_ref);
			summary_file.write("\nHet\t " + balance_het / n_het);
			summary_file.write("\nHom\t " + balance_hom / n_hom);

			//mismatch counts
			summary_file.write("\nReference Mismatch Counts\n--\n");

			for (String key : mismatch_counts.keySet()) {
				Integer h2 = mismatch_counts.get(key);
				summary_file.write(key.charAt(0) + "->" + key.charAt(1) + "\t" + h2.toString() + "\t" + ((double) h2) / mismatch_total * 100 + "\n");
			}

			summary_file.close();


			///validation table
			validation_file.print("\t");
			for (CallType c2 : CallType.values()) {
				validation_file.print(c2);
				validation_file.print("\t");
			}
			validation_file.println();

			for (CallType c1 : CallType.values()) {
				validation_file.print(c1);
				validation_file.print("\t");
				for (CallType c2 : CallType.values()) {
					SNPCallPair pair = new SNPCallPair(c1, c2);
					Long count = calltable.get(pair);
					if (count == null)
						validation_file.print(0);
					else
						validation_file.print(count);
					validation_file.print("\t");
				}
				validation_file.println();
			}
			validation_file.close();

			//strand coverage
			for (int x = 0; x < 50; x++) {
				for (int y = 0; y < 50; y++)
					strand_coverage_file.write(strand_coverage[x][y] + "\t");
				strand_coverage_file.write("\n");
			}
			strand_coverage_file.close();

			//validation breakdown by coverage
			for (int x = 0; x < Math.min(max_coverage, 65535); x++) {
				validation_coverage_file.printf("%d\t%d\t%d\t%d\n", x, calls_at_coverage[x], wrongcalls_at_coverage[x], missedcalls_at_coverage[x]);
			}

			validation_coverage_file.close();

			//validation , coverage breakdown
			for (CallType c1 : CallType.values()) {
				for (CallType c2 : CallType.values()) {
					SNPCallPair pair = new SNPCallPair(c1, c2);
					TreeMap<Integer, Long> cov = calltable_validation.get(pair);
					if (cov != null) {
						File v = new File(OUTPUT + "_summary" + File.separator + c1.toString() + c2.toString());
						PrintWriter vwriter = null;
						try {
							vwriter = new PrintWriter(v);
						} catch (FileNotFoundException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						for (Integer i : cov.keySet()) {
							vwriter.print(i + "\t");

							Long l = cov.get(i);
							if (l != null)
								vwriter.print(l);
							else
								vwriter.print(0);
							vwriter.println();
						}
						vwriter.close();
					}
				}
				validation_file.println();
			}
			validation_file.close();


		}

		return 0;
	}

	private double CalculateAlleleBalance(List<MappedBaseInfo> mappedBases,
										  Character referenceNucleotide) {
		double refcount = 0;
		double altcount = 0;

		for (MappedBaseInfo m : mappedBases) {
			if (m.nucleotide == referenceNucleotide)
				refcount += 1 - Math.pow(10, (-m.quality / 10));
			else
				altcount += 1 - Math.pow(10, (-m.quality / 10));

		}

		return (altcount / (altcount + refcount));

	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.exit(new SolSNP().instanceMain(args));
	}


	public class SNPCallPair implements Cloneable {
		public SNPCallPair(CallType one, CallType two) {
			call1 = one;
			call2 = two;
		}

		public CallType call1;
		public CallType call2;

		@Override
		public Object clone() {
			return new SNPCallPair(call1, call2);
		}

		@Override
		public int hashCode() {
			return 0;
		}


		@Override
		public boolean equals(Object obj) {
			return (call1.equals(((SNPCallPair) obj).call1) && call2.equals(((SNPCallPair) obj).call2));
		}

		@Override
		public String toString() {
			return call1.toString() + "/" + call2.toString();
		}
	}

}


