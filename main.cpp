#include <assert.h>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>

// there may be some issue here. with a different include order
// (assembler/caller/smoother at the end), compilation fails
#include "assembler.hpp"
#include "caller.hpp"
#include "chromosomes.hpp"
#include "config.hpp"
#include "lprint.hpp"
#include "ping_pong.hpp"
#include "smoother.hpp"

using namespace std;

const string version = "v1.0.5";

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void create_workdir() {
  auto c = Configuration::getInstance();
  struct stat info;
  // lprint({c->workdir}, 0);
  if (stat(c->workdir.c_str(), &info) != 0) {
    // lprint({"Working directory does not exist. creating.."}, 0);
    int check = mkdir(c->workdir.c_str(), 0777);
    if (check != 0) {
      // lprint({"Failed to create output directory, aborting.."}, 2);
      exit(EXIT_FAILURE);
    }
  }
}

int main(int argc, char **argv) {
  time_t t;
  time(&t);
  auto c = Configuration::getInstance();

  if (argc == 1) {
    c->print_help("");
    exit(EXIT_FAILURE);
  }

  c->parse(argc, argv);

  if (c->version) {
    cout << "SVDSS, " << version << endl;
    exit(EXIT_SUCCESS);
  }

  if (c->help) {
    c->print_help(argv[1]);
    exit(EXIT_SUCCESS);
  }

  // cerr << "SVDSS, Structural Variant Discovery from Sample-specific Strings."
  // << endl; cerr << "Mode: " << argv[1] << endl;
  if (strcmp(argv[1], "call") == 0) {
    if (c->reference == "" || c->bam == "" || c->workdir == "") {
      c->print_help(argv[1]);
      exit(EXIT_FAILURE);
    }
    auto caller = new Caller();
    caller->run();
  } else if (strcmp(argv[1], "index") == 0) {
    if (c->reference == "" || c->index == "") {
      c->print_help(argv[1]);
      exit(EXIT_FAILURE);
    }
    auto pingpong = new PingPong();
    pingpong->index();
  } else if (strcmp(argv[1], "search") == 0) {
    if (c->index == "" || c->bam == "" || c->workdir == "") {
      c->print_help(argv[1]);
      exit(EXIT_FAILURE);
    }
    create_workdir();
    auto pingpong = new PingPong();
    pingpong->search();
  } else if (strcmp(argv[1], "assemble") == 0) {
    auto assembler = new Assembler();
    assembler->run();
  } else if (strcmp(argv[1], "smooth") == 0) {
    if (c->reference == "" || c->bam == "") {
      c->print_help(argv[1]);
      exit(EXIT_FAILURE);
    }
    auto smoother = new Smoother();
    smoother->run();
  } else {
    c->print_help("");
    exit(EXIT_FAILURE);
  }
  time_t s;
  time(&s);
  lprint({"Complete. Runtime:", to_string(s - t) + " seconds."});

  return 0;
}
