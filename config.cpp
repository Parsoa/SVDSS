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
  // clang-format off
  parser.add_options()
    ("bam", "", cxxopts::value<std::string>())
    ("sfs", "", cxxopts::value<std::string>())
    ("poa", "", cxxopts::value<std::string>())
    ("clusters", "", cxxopts::value<std::string>())
    ("index", "", cxxopts::value<std::string>())
    ("fastx", "", cxxopts::value<std::string>())
    ("reference", "", cxxopts::value<std::string>())
    ("append", "", cxxopts::value<std::string>())
    ("threads", "", cxxopts::value<int>())
    ("bsize", "", cxxopts::value<int>())
    ("omax", "", cxxopts::value<int>())
    ("min-sv-length", "", cxxopts::value<int>())
    ("min-mapq", "", cxxopts::value<int>())
    ("min-cluster-weight", "", cxxopts::value<int>())
    ("clipped", "", cxxopts::value<bool>()->default_value("false"))
    ("noref", "", cxxopts::value<bool>()->default_value("false"))
    ("noassemble", "", cxxopts::value<bool>()->default_value("false"))
    ("noputative", "", cxxopts::value<bool>()->default_value("false"))
    ("binary", "", cxxopts::value<bool>()->default_value("false"))
    ("version", "Print version information.")
    ("h,help", "Print this help.")
    ("l", "", cxxopts::value<float>())
    ("acc", "", cxxopts::value<float>())
    ("verbose", "", cxxopts::value<bool>()->default_value("false"));
  // clang-format on
}

void Configuration::parse(int argc, char **argv) {
  auto results = parser.parse(argc, argv);
  if (results.count("bam"))
    bam = results["bam"].as<std::string>();
  if (results.count("sfs"))
    sfs = results["sfs"].as<std::string>();
  if (results.count("poa"))
    poa = results["poa"].as<std::string>();
  if (results.count("clusters"))
    clusters = results["clusters"].as<std::string>();
  if (results.count("index"))
    index = results["index"].as<std::string>();
  if (results.count("fastx"))
    fastq = results["fastx"].as<std::string>(); // FIXME: use fastx here and in pingpong
  if (results.count("overlap"))
    overlap = results["overlap"].as<int>();
  if (results.count("bsize"))
    batch_size = results["bsize"].as<int>();
  if (results.count("omax"))
    max_output = results["omax"].as<int>();
  if (results.count("threads"))
    threads = results["threads"].as<int>();
  if (results.count("append"))
    append = results["append"].as<std::string>();
  if (results.count("reference"))
    reference = results["reference"].as<std::string>();
  if (results.count("min-sv-length"))
    min_sv_length = max(25, results["min-sv-length"].as<int>());
  if (results.count("min-cluster-weight"))
    min_cluster_weight = results["min-cluster-weight"].as<int>();
  if (results.count("min-mapq"))
    min_mapq = results["min-mapq"].as<int>();
  if (results.count("l"))
    min_ratio = results["l"].as<float>();
  if (results.count("acc"))
    al_accuracy = results["acc"].as<float>();
  binary = results["binary"].as<bool>();
  clipped = results["clipped"].as<bool>();
  noref = results["noref"].as<bool>();
  assemble = !(results["noassemble"].as<bool>());
  putative = !(results["noputative"].as<bool>());
  version = results["version"].as<bool>();
  verbose = results["verbose"].as<bool>();
  help = results["help"].as<bool>();

  batch_size = (batch_size / threads) * threads;
}
