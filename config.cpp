#include "config.hpp"

Configuration *Configuration::instance = nullptr;

Configuration *Configuration::getInstance() {
  if (instance == nullptr) {
    instance = new Configuration();
  }
  return instance;
}

void Configuration::print_help(const string &mode) const {
  if (mode.compare("index") == 0) {
    cerr << INDEX_USAGE_MESSAGE << endl;
  } else if (mode.compare("smooth") == 0) {
    cerr << SMOOTH_USAGE_MESSAGE << endl;
  } else if (mode.compare("search") == 0) {
    cerr << SEARCH_USAGE_MESSAGE << endl;
  } else if (mode.compare("call") == 0) {
    cerr << CALL_USAGE_MESSAGE << endl;
  } else {
    cerr << MAIN_USAGE_MESSAGE << endl;
  }
}

Configuration::Configuration()
    : parser(
          "SVDSS, Structural Variant Discovery from Sample-specific Strings") {
  parser.add_options()("bed", "", cxxopts::value<std::string>())(
      "bam", "", cxxopts::value<std::string>())("vcf", "",
                                                cxxopts::value<std::string>())(
      "type", "", cxxopts::value<std::string>())("index", "",
                                                 cxxopts::value<std::string>())(
      "fastq", "", cxxopts::value<std::string>())(
      "target", "", cxxopts::value<std::string>())(
      "sfsbam", "", cxxopts::value<std::string>())(
      "append", "", cxxopts::value<std::string>())(
      "workdir", "", cxxopts::value<std::string>())(
      "reference", "", cxxopts::value<std::string>())("cutoff", "",
                                                      cxxopts::value<int>())(
      "coverage", "", cxxopts::value<int>())("t,threads", "",
                                             cxxopts::value<int>())(
      "min-sv-length", "", cxxopts::value<int>())("min-cluster-weight", "",
                                                  cxxopts::value<int>())(
      "clipped", "", cxxopts::value<bool>()->default_value("false"))(
      "noassemble", "", cxxopts::value<bool>()->default_value("false"))(
      "noputative", "", cxxopts::value<bool>()->default_value("false"))(
      "b,binary", "", cxxopts::value<bool>()->default_value("false"))(
      "aggregate", "", cxxopts::value<bool>()->default_value("false"))(
      "selective", "", cxxopts::value<bool>()->default_value("true"))(
      "version", "Print version information.")("h,help", "Print this help.")(
      "l", "", cxxopts::value<float>())("acc", "", cxxopts::value<float>())(
      "verbose", "", cxxopts::value<bool>()->default_value("false"));
}

void Configuration::parse(int argc, char **argv) {
  auto results = parser.parse(argc, argv);
  if (results.count("vcf")) {
    vcf = results["vcf"].as<std::string>();
  }
  bam = "";
  if (results.count("bam")) {
    bam = results["bam"].as<std::string>();
  }
  if (results.count("sfsbam")) {
    sfsbam = results["sfsbam"].as<std::string>();
  }
  if (results.count("bed")) {
    bed = results["bed"].as<std::string>();
  }
  if (results.count("type")) {
    type = results["type"].as<std::string>();
  }
  if (results.count("index")) {
    index = results["index"].as<std::string>();
  }
  fastq = "";
  if (results.count("fastq")) {
    fastq = results["fastq"].as<std::string>();
  }
  target = "";
  if (results.count("target")) {
    target = results["target"].as<std::string>();
  }
  if (results.count("cutoff")) {
    cutoff = results["cutoff"].as<int>();
  }
  if (results.count("overlap")) {
    overlap = results["overlap"].as<int>();
  }
  if (results.count("threads")) {
    threads = results["threads"].as<int>();
  }
  if (results.count("workdir")) {
    workdir = results["workdir"].as<std::string>();
  }
  if (results.count("append")) {
    append = results["append"].as<std::string>();
  } else {
    append = "";
  }
  if (results.count("coverage")) {
    coverage = results["coverage"].as<int>();
  }
  if (results.count("reference")) {
    reference = results["reference"].as<std::string>();
  }
  if (results.count("min-sv-length")) {
    min_sv_length = max(25, results["min-sv-length"].as<int>());
  }
  if (results.count("min-cluster-weight")) {
    min_cluster_weight = results["min-cluster-weight"].as<int>();
  }
  if (results.count("l")) {
    min_ratio = results["l"].as<float>();
  }
  if (results.count("acc")) {
    al_accuracy = results["acc"].as<float>();
  }
  binary = results["binary"].as<bool>();
  clipped = results["clipped"].as<bool>();
  assemble = !(results["noassemble"].as<bool>());
  putative = !(results["noputative"].as<bool>());
  aggregate = results["aggregate"].as<bool>();
  selective = results["selective"].as<bool>();
  version = results["version"].as<bool>();
  verbose = results["verbose"].as<bool>();
  help = results["help"].as<bool>();
}
