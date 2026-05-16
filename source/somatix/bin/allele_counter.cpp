// allele_counter.cpp
// Compile:
//   g++ -O3 -std=c++17 allele_counter.cpp -o allele_counter -lhts -lz -lbz2 -llzma -lcurl -lpthread
//   g++ -O3 -std=c++17 allele_counter.cpp -o allele_counter -I /home/colin/biotools/htslib-1.23.1  -L/home/colin/biotools/htslib-1.23.1   -lhts -lz -lbz2 -llzma -lcurl -lpthread
#include <htslib/sam.h>
#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Args {
    std::string bam_path;
    std::string ref_path;
    std::string region;

    int min_mapq = 20;
    int min_baseq = 10;
    int min_alt = 2;

    // total_coverage = A+C+G+T+INS+DEL
    int min_total_coverage = 10;

    // vaf = alt_count / total_coverage
    double min_vaf = 0.05;

    // cap using A+C+G+T+N+INS+DEL
    int max_depth = 5000;

    int excl_flags = 2316;
};

struct Count {
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
    int N = 0;
    int INS = 0;
    int DEL = 0;
};

static void print_usage() {
    std::cerr
        << "Usage:\n"
        << "  allele_counter --bam sample.bam --ref ref.fa --region chr1:1-1000000 "
        << "[--min-mapq 20] [--min-baseq 10] [--min-alt 2] "
        << "[--min-total-coverage 10] [--min-vaf 0.05] "
        << "[--max-depth 5000] [--excl-flags 2316]\n";
}

static Args parse_args(int argc, char** argv) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];

        auto need_value = [&](const std::string& k) {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + k);
            }
            return std::string(argv[++i]);
        };

        if (key == "--bam") {
            args.bam_path = need_value(key);
        } else if (key == "--ref") {
            args.ref_path = need_value(key);
        } else if (key == "--region") {
            args.region = need_value(key);
        } else if (key == "--min-mapq" || key == "--min-MQ") {
            args.min_mapq = std::stoi(need_value(key));
        } else if (key == "--min-baseq" || key == "--min-BQ") {
            args.min_baseq = std::stoi(need_value(key));
        } else if (key == "--min-alt") {
            args.min_alt = std::stoi(need_value(key));
        } else if (key == "--min-total-coverage") {
            args.min_total_coverage = std::stoi(need_value(key));
        } else if (key == "--min-vaf") {
            args.min_vaf = std::stod(need_value(key));
        } else if (key == "--max-depth") {
            args.max_depth = std::stoi(need_value(key));
        } else if (key == "--excl-flags") {
            args.excl_flags = std::stoi(need_value(key));
        } else if (key == "--help" || key == "-h") {
            print_usage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }

    if (args.bam_path.empty() || args.ref_path.empty() || args.region.empty()) {
        print_usage();
        throw std::runtime_error("Missing required arguments.");
    }

    return args;
}

static char decode_base(uint8_t code) {
    static const char* seq_nt16_str = "=ACMGRSVTWYHKDBN";
    char b = seq_nt16_str[code & 0xF];

    switch (b) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return b;
        default:
            return 'N';
    }
}

static void add_base_count(Count& c, char base) {
    switch (base) {
        case 'A': c.A++; break;
        case 'C': c.C++; break;
        case 'G': c.G++; break;
        case 'T': c.T++; break;
        default:  c.N++; break;
    }
}

static int get_base_count(const Count& c, char base) {
    switch (base) {
        case 'A': return c.A;
        case 'C': return c.C;
        case 'G': return c.G;
        case 'T': return c.T;
        default:  return 0;
    }
}

static int get_total_depth_all_with_N(const Count& c) {
    return c.A + c.C + c.G + c.T + c.N + c.INS + c.DEL;
}

static int get_total_coverage_6types(const Count& c) {
    return c.A + c.C + c.G + c.T + c.INS + c.DEL;
}

static std::string best_alt(const Count& c, char ref, int& alt_count) {
    // Select ALT by maximum read count.
    // If tied, use priority:
    //  A > C > G > T > INS > DEL >

    alt_count = 0;
    std::string alt = "N";

    auto consider = [&](const std::string& allele, int count) {
        if (count > alt_count) {
            alt_count = count;
            alt = allele;
        }
    };

    // Tie priority is controlled by this order.
    // Since we only update when count > alt_count, the first allele
    // with the maximum count is kept.
    

    if (ref != 'A') consider("A", c.A);
    if (ref != 'T') consider("T", c.T);
    if (ref != 'C') consider("C", c.C);
    if (ref != 'G') consider("G", c.G);
    

    consider("INS", c.INS);
    consider("DEL", c.DEL);
    
    return alt;
}


int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);

        samFile* in = sam_open(args.bam_path.c_str(), "r");
        if (!in) {
            throw std::runtime_error("Cannot open BAM/CRAM: " + args.bam_path);
        }

        bam_hdr_t* hdr = sam_hdr_read(in);
        if (!hdr) {
            sam_close(in);
            throw std::runtime_error("Cannot read BAM header.");
        }

        hts_idx_t* idx = sam_index_load(in, args.bam_path.c_str());
        if (!idx) {
            bam_hdr_destroy(hdr);
            sam_close(in);
            throw std::runtime_error("Cannot load BAM index. Please run samtools index first.");
        }

        faidx_t* fai = fai_load(args.ref_path.c_str());
        if (!fai) {
            hts_idx_destroy(idx);
            bam_hdr_destroy(hdr);
            sam_close(in);
            throw std::runtime_error("Cannot load FASTA index. Please run samtools faidx first.");
        }

        hts_itr_t* itr = sam_itr_querys(idx, hdr, args.region.c_str());
        if (!itr) {
            fai_destroy(fai);
            hts_idx_destroy(idx);
            bam_hdr_destroy(hdr);
            sam_close(in);
            throw std::runtime_error("Invalid region or region not found: " + args.region);
        }

        int tid = itr->tid;
        int region_beg = itr->beg;
        int region_end = itr->end;

        std::string chrom = hdr->target_name[tid];
        int region_len = region_end - region_beg;

        std::vector<Count> counts(region_len);

        int ref_len = 0;
        char* ref_seq_raw = faidx_fetch_seq(
            fai,
            chrom.c_str(),
            region_beg,
            region_end - 1,
            &ref_len
        );

        if (!ref_seq_raw || ref_len <= 0) {
            throw std::runtime_error("Cannot fetch reference sequence for region.");
        }

        std::string ref_seq(ref_seq_raw, ref_len);
        free(ref_seq_raw);

        for (char& b : ref_seq) {
            b = std::toupper(static_cast<unsigned char>(b));
            if (b != 'A' && b != 'C' && b != 'G' && b != 'T') {
                b = 'N';
            }
        }

        bam1_t* read = bam_init1();

        while (sam_itr_next(in, itr, read) >= 0) {
            const bam1_core_t& core = read->core;

            if (core.flag & args.excl_flags) continue;
            if (core.qual < args.min_mapq) continue;

            uint32_t* cigar = bam_get_cigar(read);
            uint8_t* seq = bam_get_seq(read);
            uint8_t* qual = bam_get_qual(read);

            int64_t ref_pos = core.pos;
            int64_t query_pos = 0;

            for (uint32_t i = 0; i < core.n_cigar; ++i) {
                int op = bam_cigar_op(cigar[i]);
                int len = bam_cigar_oplen(cigar[i]);

                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    for (int j = 0; j < len; ++j) {
                        int64_t gpos = ref_pos + j;
                        int64_t qpos = query_pos + j;

                        if (gpos < region_beg || gpos >= region_end) {
                            continue;
                        }

                        Count& c = counts[gpos - region_beg];

                        if (get_total_depth_all_with_N(c) >= args.max_depth) {
                            continue;
                        }

                        if (qual && qual[qpos] < args.min_baseq) {
                            continue;
                        }

                        char base = decode_base(bam_seqi(seq, qpos));
                        add_base_count(c, base);
                    }

                    ref_pos += len;
                    query_pos += len;

                } else if (op == BAM_CINS) {
                    int64_t anchor_pos = ref_pos - 1;

                    if (anchor_pos >= region_beg && anchor_pos < region_end) {
                        Count& c = counts[anchor_pos - region_beg];

                        if (get_total_depth_all_with_N(c) < args.max_depth) {
                            bool pass_bq = true;

                            if (qual) {
                                pass_bq = false;
                                for (int j = 0; j < len; ++j) {
                                    if (qual[query_pos + j] >= args.min_baseq) {
                                        pass_bq = true;
                                        break;
                                    }
                                }
                            }

                            if (pass_bq) {
                                c.INS++;
                            }
                        }
                    }

                    query_pos += len;

                } else if (op == BAM_CDEL) {
                    for (int j = 0; j < len; ++j) {
                        int64_t gpos = ref_pos + j;

                        if (gpos >= region_beg && gpos < region_end) {
                            Count& c = counts[gpos - region_beg];

                            if (get_total_depth_all_with_N(c) < args.max_depth) {
                                c.DEL++;
                            }
                        }
                    }

                    ref_pos += len;

                } else if (op == BAM_CREF_SKIP) {
                    ref_pos += len;

                } else if (op == BAM_CSOFT_CLIP) {
                    query_pos += len;

                } else if (op == BAM_CHARD_CLIP) {
                    // consumes neither reference nor query
                } else if (op == BAM_CPAD) {
                    // consumes neither reference nor query
                }
            }
        }

        bam_destroy1(read);

        std::cout
            << "chrom\tpos\tref\talt\tref_count\talt_count"
            << "\ttotal_coverage\tvaf"
            << "\tA\tC\tG\tT\tINS\tDEL\n";

        for (int i = 0; i < region_len && i < ref_len; ++i) {
            const Count& c = counts[i];
            char ref = ref_seq[i];

            int alt_count = 0;
            std::string alt = best_alt(c, ref, alt_count);

            int ref_count = get_base_count(c, ref);
            int total_coverage = get_total_coverage_6types(c);

            if (total_coverage < args.min_total_coverage) {
                continue;
            }

            if (alt_count < args.min_alt) {
                continue;
            }

            double vaf = total_coverage > 0
                ? static_cast<double>(alt_count) / static_cast<double>(total_coverage)
                : 0.0;

            if (vaf < args.min_vaf) {
                continue;
            }

            std::cout
                << chrom << '\t'
                << (region_beg + i + 1) << '\t'
                << ref << '\t'
                << alt << '\t'
                << ref_count << '\t'
                << alt_count << '\t'
                << total_coverage << '\t'
                << vaf << '\t'
                << c.A << '\t'
                << c.C << '\t'
                << c.G << '\t'
                << c.T << '\t'
                << c.INS << '\t'
                << c.DEL << '\n';
        }

        hts_itr_destroy(itr);
        fai_destroy(fai);
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
