#include <ctime>

#include "caller.hpp"

using namespace std;

void Caller::run() {
  config = Configuration::getInstance();
  lprint({"PingPong SV Caller running on", to_string(config->threads),
          "threads.."});
  // load reference genome and SFS
  load_chromosomes(config->reference);
  lprint({"Loaded all chromosomes."});
  load_input_sfs();
  // call SVs from extended SFS
  vector<SV> svs;
  Extender extender = Extender(&SFSs);
  extender.run(config->threads);
  svs.insert(svs.begin(), extender.svs.begin(), extender.svs.end());
  lprint({"Predicted", to_string(svs.size()), "SVs from extended SFS."});
  interval_tree_t<int> vartree;
  for (const auto &sv : svs) {
    vartree.insert({sv.s - 1000, sv.e + 1000});
  }
  std::sort(svs.begin(), svs.end());
  // output POA alignments SAM
  string poa_path = config->workdir + "/poa.sam";
  lprint({"Outputting POA alignments to", poa_path + ".."});
  osam.open(poa_path);
  osam << "@HD\tVN:1.4" << endl;
  for (int i = 0; i < chromosomes.size(); ++i) {
    osam << "@SQ\tSN:" << chromosomes[i] << "\t"
         << "LN:" << strlen(chromosome_seqs[chromosomes[i]]) << endl;
  }
  for (int j = 0; j < extender.alignments.size(); j++) {
    const auto &c = extender.alignments[j];
    osam << c.chrom << ":" << c.s + 1 << "-" << c.e + 1 << "\t"
         << "0"
         << "\t" << c.chrom << "\t" << c.s + 1 << "\t"
         << "60"
         << "\t" << c.cigar << "\t"
         << "*"
         << "\t"
         << "0"
         << "\t"
         << "0"
         << "\t" << c.seq << "\t"
         << "*" << endl;
  }
  osam.close();
  // output SV calls
  string vcf_path = config->workdir + "/svs_poa.vcf";
  lprint({"Exporting", to_string(svs.size()), "SV calls to", vcf_path + ".."});
  ovcf.open(vcf_path);
  print_vcf_header();
  for (const SV &sv : svs) {
    ovcf << sv << endl;
  }
  ovcf.close();
  //
  if (config->clipped) {
    vector<SV> clipped_svs;
    Clipper clipper(extender.clips);
    clipper.call(config->threads, vartree);
    int s = 0;
    for (int i = 0; i < config->threads; i++) {
      s += clipper._p_svs[i].size();
      clipped_svs.insert(svs.begin(), clipper._p_svs[i].begin(),
                         clipper._p_svs[i].end());
    }
    lprint({"Predicted", to_string(s), "SVs from clipped SFS."});
    string vcf_path = config->workdir + "/svs_clipped.vcf";
    lprint({"Exporting", to_string(clipped_svs.size()), "SV calls to",
            vcf_path + ".."});
    ovcf.open(vcf_path);
    print_vcf_header();
    for (const SV &sv : clipped_svs) {
      ovcf << sv << endl;
    }
    ovcf.close();
  }
}

void Caller::load_input_sfs() {
  int threads = config->threads;
  int num_batches = 0;
  if (auto dir = opendir(config->workdir.c_str())) {
    while (auto f = readdir(dir)) {
      string fn = f->d_name; // FIXME
      if (!f->d_name || f->d_name[0] == '.')
        continue; // Skip everything that starts with a dot
      if (fn.substr(0, 5).compare("solut") == 0)
        ++num_batches;
    }
    closedir(dir);
  }

  vector<unordered_map<string, vector<SFS>>> _SFSs(num_batches);
  int num_threads = num_batches < threads ? num_batches : threads;
  lprint({"Loading assmbled SFS.."});
#pragma omp parallel for num_threads(num_threads)
  for (int j = 0; j < num_batches; j++) {
    string s_j = std::to_string(j);
    string inpath =
        config->workdir + "/solution_batch_" + s_j + ".assembled.sfs";
    // cout << "[I] Loading SFS from " << inpath << endl ;
    ifstream inf(inpath);
    string line;
    if (inf.is_open()) {
      string info[4];
      string read_name;
      while (getline(inf, line)) {
        stringstream ssin(line);
        int i = 0;
        while (ssin.good() && i < 4) {
          ssin >> info[i++];
        }
        if (info[0].compare("*") != 0) {
          read_name = info[0];
          _SFSs[j][read_name] = vector<SFS>();
        }
        _SFSs[j][read_name].push_back(
            SFS(stoi(info[1]), stoi(info[2]), stoi(info[3]), true));
      }
    }
  }
  int r = 0;
  int c = 0;
  for (int j = 0; j < num_batches; j++) {
    // lprint({"Batch", to_string(j), "with", to_string(_SFSs[j].size()),
    // "strings."});
    r += _SFSs[j].size();
    SFSs.insert(_SFSs[j].begin(), _SFSs[j].end());
    for (auto &read : _SFSs[j]) {
      c += read.second.size();
    }
  }
  lprint({"Loaded", to_string(c), "SFS strings on", to_string(r), "reads."});
}

void Caller::print_vcf_header() {
  ovcf << "##fileformat=VCFv4.2" << endl;
  ovcf << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
          "data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/"
          "hg38.no_alt.fa.gz"
       << endl;
  for (int i = 0; i < chromosomes.size(); ++i) {
    ovcf << "##contig=<ID=" << chromosomes[i]
         << ",length=" << strlen(chromosome_seqs[chromosomes[i]]) << ">"
         << endl;
  }
  ovcf << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
  ovcf << "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant "
          "class\">"
       << endl;
  ovcf << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">"
       << endl;
  ovcf << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in "
          "length between REF and ALT alleles\">"
       << endl;
  ovcf << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of "
          "the variant described in this record\">"
       << endl;
  ovcf << "##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number of "
          "alignments supporting this record\">"
       << endl;
  ovcf << "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total number of "
          "alignments covering this locus\">"
       << endl;
  ovcf << "##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Alignment score\">"
       << endl;
  ovcf << "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of "
          "variations on same consensus\">"
       << endl;
  ovcf << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise "
          "structural variation\">"
       << endl;
  ovcf << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR of "
          "consensus\">"
       << endl;
  ovcf << "##INFO=<ID=READS,Number=.,Type=String,Description=\"Reads "
          "identifiers supporting the call\">"
       << endl;
  ovcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
       << endl;
  ovcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEFAULT"
       << endl;
}
