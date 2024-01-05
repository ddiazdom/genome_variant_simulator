#include <iostream>
#include <zlib.h>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include "kseq.h"
#include "CLI11.hpp"

KSEQ_INIT(int, read)

struct arguments{
    std::string input_file;
    std::string output_file;
    size_t n_genomes;
    float var_frac=0.01;
    size_t max_edit_len=10;
};

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option * opt) const override {
        std::stringstream out;
        if(!opt->get_option_text().empty()) {
            out << " " << opt->get_option_text();
        } else {
            if(opt->get_type_size() != 0) {
                if(!opt->get_type_name().empty()) out << " " << get_label(opt->get_type_name());
                //if(!opt->get_default_str().empty()) out << "=" << opt->get_default_str();
                if(opt->get_expected_max() == CLI::detail::expected_max_vector_size) out << " ...";
                else if(opt->get_expected_min() > 1) out << " x " << opt->get_expected();

                if(opt->get_required())
                    out << " " << get_label("REQUIRED");
            }
        }
        return out.str();
    }
};

int parse_app(CLI::App& app, struct arguments& args) {

    auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(29);
    app.formatter(fmt);

    app.add_option("INPUT_GENOME", args.input_file, "input file FASTX with the genome")->check(CLI::ExistingFile)->required();
    app.add_option("N_GENOMES", args.n_genomes, "number of genomes to simulate")->check(CLI::PositiveNumber)->required();
    app.add_option("OUTPUT_FILE", args.output_file, "output file where the simulated genomes will be stored");
    app.add_option("-f,--var-frac", args.var_frac, "fraction of the genome size to be modified (def. 0.01)")->check(CLI::Range(0.0, 0.9999));
    app.add_option("-d,--max-edit-len", args.max_edit_len, "maximum length of a variant (def. 10)");

    return 0;
}

void create_variants(std::vector<std::string>& genome, size_t g_len, float var_frac,
                     long max_edit_len, size_t n_genomes,
                     std::string& output_file){

    char dna_symbols[8] = {'a', 'c', 'g', 't', 'A', 'C', 'T', 'G'};
    std::string buffer;
    buffer.resize(max_edit_len);
    std::ofstream ofs(output_file, std::ios::binary | std::ios::out);

    for(size_t i=0;i<n_genomes;i++) {
        //event 0 = SNP
        //event 1 = insertion
        //event 2 = deletion
        //event 3 = inversion
        //event 4 = translocation
        std::random_device dev;
        std::mt19937 rng(dev());
        std::discrete_distribution<int> distribution {90,3,3,2,2};
        std::vector<std::string> new_genome = genome;
        size_t events[5] = {0};
        size_t sym_event[5] = {0};

        long n_vars = ceil(float(g_len)*var_frac);
        while (n_vars > 0) {
            int var_event = distribution(rng);
            events[var_event]++;

            if (var_event == 0) {//SNP
                n_vars--;
                size_t chrom = rand() % new_genome.size();
                size_t pos = rand() % new_genome[chrom].size();
                char new_sym = dna_symbols[rand() % 8];
                while (new_sym == new_genome[chrom][pos]) {
                    new_sym = dna_symbols[rand() % 8];
                }
                new_genome[chrom][pos] = new_sym;
                sym_event[var_event]++;
            } else {
                size_t chrom = rand() % new_genome.size();
                long ed_len = 1+rand() % long(max_edit_len);

                assert(ed_len <= max_edit_len);
                long pos = rand() % (long(new_genome[chrom].size()) - ed_len);
                assert(pos + ed_len <= long(new_genome[chrom].size()));

                if (var_event == 1) {//insertion
                    for(long j = 0; j < ed_len; j++) {//random string
                        buffer[j] = dna_symbols[rand() % 8];
                    }
                    new_genome[chrom].insert(pos, buffer, 0, ed_len);
                    n_vars -= ed_len;
                    sym_event[var_event]+=ed_len;
                } else if (var_event == 2) {//deletion
                    new_genome[chrom].erase(pos, ed_len);
                    n_vars -= ed_len;
                    sym_event[var_event]+=ed_len;
                } else if (var_event == 3) {//inversion
                    assert(pos + ed_len <= long(new_genome[chrom].size()));
                    std::reverse(new_genome[chrom].begin() + pos, new_genome[chrom].begin() + pos + ed_len);
                    n_vars -= ed_len;
                    sym_event[var_event]+=ed_len;
                } else if (var_event == 4) {//translocation

                    size_t chrom_2 = rand() % new_genome.size();
                    long pos_2 = rand() % (long(new_genome[chrom_2].size()) - ed_len);
                    while(chrom_2 == chrom && pos == pos_2) {
                        chrom_2 = rand() % new_genome.size();
                        pos_2 = rand() % (long(new_genome[chrom_2].size()) - ed_len);
                    }
                    assert(pos_2 + ed_len <= long(new_genome[chrom_2].size()));
                    char *ptr1 = new_genome[chrom].data() + pos;
                    char *ptr2 = new_genome[chrom_2].data() + pos_2;

                    memcpy(buffer.data(), ptr1, ed_len);
                    memcpy(ptr1, ptr2, ed_len);
                    memcpy(ptr2, buffer.data(), ed_len);
                    n_vars -= 2 * ed_len;
                    sym_event[var_event]+=2*ed_len;
                } else {
                    std::cout << "Error" << std::endl;
                }
            }
        }
        size_t n_mods=0;
        for(unsigned long j : sym_event){
            n_mods+=j;
        }

        std::cout<<"Genome "<<i+1<<" stats :"<<std::endl;
        std::cout << "  Modified symbols: " << n_mods << std::endl;
        std::cout << "  SNPs " << events[0] << ", insertions " << events[1] << ", deletions " << events[2] << ", inversions "
                  << events[3] << ", translocation " << events[4] << std::endl;

        for(auto & chrom : new_genome){
            chrom.push_back('\n');
            ofs.write(chrom.data(), long(chrom.size()));
        }
    }
    ofs.close();
}

int main(int argc, char ** argv) {

    CLI::App app("Genome variant simulator");
    arguments argu;
    parse_app(app, argu);
    CLI11_PARSE(app, argc, argv);
    FILE* fp;
    kseq_t *seq;

    fp = fopen(argu.input_file.c_str(), "r");
    seq = kseq_init(fileno(fp));

    std::vector<std::string> genome;
    size_t g_len=0;
    while (kseq_read(seq) >= 0){
        genome.emplace_back(seq->seq.s);
        g_len +=seq->seq.l;
    }
    genome.shrink_to_fit();
    kseq_destroy(seq);
    fclose(fp);

    create_variants(genome, g_len, argu.var_frac, long(argu.max_edit_len), argu.n_genomes, argu.output_file);
    return 0;
}
