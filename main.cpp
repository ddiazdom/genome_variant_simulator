#include <iostream>
#include <zlib.h>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <cassert>
#include "kseq.h"
#include "CLI11.hpp"

KSEQ_INIT(int, read)

char dna_symbols[8] = {'a', 'c', 'g', 't', 'A', 'C', 'T', 'G'};

struct arguments{
    std::string input_file;
    std::string output_file;
    size_t n_genomes=0;
    float var_frac=0.01;
    size_t max_edit_len=10;
};

struct variant{
    int type;
    long pos_1{};
    long pos_2{};
    size_t chrom_1{};
    size_t chrom_2{};
    long len{};
    char alt_allele{};
    std::string alt_pat;
    std::vector<bool> occ;
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

void apply_variants(std::vector<std::string>& genome, size_t n_genomes, std::string& output_file,
                    size_t max_edit_len, std::vector<variant>& variants){
    std::string buffer;
    buffer.resize(max_edit_len);
    std::ofstream ofs(output_file, std::ios::binary | std::ios::out);
    char *ptr1, *ptr2;
    std::vector<long> chrom_offset(genome.size(), 0);

    for(size_t i=0;i<n_genomes;i++) {
        std::vector<std::string> new_genome = genome;

        for(auto & variant : variants){
            if(variant.occ[i]){
                switch (variant.type) {
                    case 0://SNP
                        new_genome[variant.chrom_1][variant.pos_1] = variant.alt_allele;
                        break;
                    case 1://translocation
                        ptr1 = new_genome[variant.chrom_1].data() + variant.pos_1;
                        ptr2 = new_genome[variant.chrom_2].data() + variant.pos_2;
                        memcpy(buffer.data(), ptr1, variant.len);
                        memcpy(ptr1, ptr2, variant.len);
                        memcpy(ptr2, buffer.data(), variant.len);
                        break;
                    case 2://inversion
                        std::reverse(new_genome[variant.chrom_1].begin() + variant.pos_1,
                                     new_genome[variant.chrom_1].begin() + variant.pos_1 + variant.len);
                        break;
                    case 3://deletion
                        //TODO I have to fix this because it is not working
                        new_genome[variant.chrom_1].erase(variant.pos_1, variant.len);
                        break;
                    case 4://insertion
                        //TODO I have to fix this because it is not working
                        new_genome[variant.chrom_1].insert(variant.pos_1, variant.alt_pat, 0, variant.len);
                        break;
                    default:
                        std::cout<<"Error"<<std::endl;
                        break;
                }
            }
        }

        for(auto & chrom : new_genome){
            chrom.push_back('\n');
            ofs.write(chrom.data(), long(chrom.size()));
        }
    }
}

void create_variants(std::vector<std::string>& genome, size_t g_len, float var_frac, long max_edit_len,
                     size_t n_genomes, std::string& output_file){

    long n_vars = ceil(float(g_len)*var_frac);

    //event 0 = SNP
    //event 1 = insertion
    //event 2 = deletion
    //event 3 = inversion
    //event 4 = translocation
    std::random_device dev;
    std::mt19937 rng(dev());
    std::discrete_distribution<int> distribution {94,3,2,1,0};
    std::vector<std::string> new_genome = genome;

    size_t events[5] = {0};
    size_t sym_event[5] = {0};

    std::vector<variant> genome_variants;
    while (n_vars > 0) {

        variant var;
        var.type = distribution(rng);
        events[var.type]++;

        if (var.type == 0){//SNP
            n_vars--;
            var.chrom_1 = rand() % new_genome.size();
            var.pos_1 = rand() % new_genome[var.chrom_1].size();

            var.alt_allele = dna_symbols[rand() % 8];
            while (var.alt_allele == new_genome[var.chrom_1][var.pos_1]) {
                var.alt_allele = dna_symbols[rand() % 8];
            }
            sym_event[var.type]++;
        } else {

            var.chrom_1 = rand() % new_genome.size();
            var.len = 1+rand() % long(max_edit_len);
            assert(var.len <= max_edit_len);

            var.pos_1 = rand() % (long(new_genome[var.chrom_1].size()) - var.len);
            assert(var.pos_1 + var.len <= long(new_genome[var.chrom_1].size()));

            if (var.type == 1) {//translocation
                var.chrom_2 = rand() % new_genome.size();
                var.pos_2 = rand() % (long(new_genome[var.chrom_2].size()) - var.len);

                while(var.chrom_2 == var.chrom_1 && var.pos_1 == var.pos_2) {
                    var.chrom_2 = rand() % new_genome.size();
                    var.pos_2 = rand() % (long(new_genome[var.chrom_2].size()) - var.len);
                }
                assert(var.pos_2 + var.len <= long(new_genome[var.chrom_2].size()));

                n_vars -= 2 * var.len;
                sym_event[var.type]+=2*var.len;
            } else if (var.type==2 || var.type == 3) {//inversion and deletion
                n_vars -= var.len;
                sym_event[var.type]+=var.len;
            } else if (var.type == 4) {//insertion
                var.alt_pat.resize(var.len);
                for(long j = 0; j < var.len; j++) {//random string
                    var.alt_pat[j] = dna_symbols[rand() % 8];
                }
                n_vars -= var.len;
                sym_event[var.type]+=var.len;
            } else {
                std::cout << "Error" << std::endl;
            }
        }

        double ref = rand() % 100;
        double alt = 100-ref;
        std::random_device dev_2;
        std::mt19937 rng2(dev_2());
        std::discrete_distribution<int> all_freq{ref, alt};
        var.occ.resize(n_genomes);
        for(size_t j=0;j<n_genomes;j++){
            var.occ[j] = all_freq(rng2);
        }
        genome_variants.emplace_back(std::move(var));
    }
    std::cout <<"SNPs " << events[0] <<", translocation " << events[1] << ", inversions " << events[2] << ", deletions "
              << events[3] << ", insertions " << events[4] << std::endl;

    std::sort(genome_variants.begin(), genome_variants.end(), [](const variant& a, const variant& b){
        if(a.type!=b.type){
            return a.type<b.type;
        }else{
            if(a.chrom_1!=b.chrom_1){
                return a.chrom_1<b.chrom_1;
            }else{
                return a.pos_1>b.pos_1;
            }
        }
    });

    /*for(size_t j=0;j<genome_variants.size();j++){
        int counts =0;
        for(size_t k=0;k<genome_variants[j].occ.size();k++){
            counts+=genome_variants[j].occ[k];
        }
        std::cout<<genome_variants[j].type<<" : "<<genome_variants[j].chrom_1<<" "<<genome_variants[j].pos_1<<" -> "<<float(counts)/float(n_genomes)<<std::endl;
    }*/
    std::cout<<" Creating the genomes "<<std::endl;
    apply_variants(genome, n_genomes, output_file, max_edit_len, genome_variants);
    std::cout<<" Genomes stored in "<<output_file<<std::endl;
}



void insert_variants(std::vector<std::string>& genome, size_t g_len, float var_frac,
                     long max_edit_len, size_t n_genomes,
                     std::string& output_file){

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
