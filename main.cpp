#include <ctime>
#include <string>
#include <sys/stat.h>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "caller.hpp"
#include "config.hpp"
#include "ping_pong.hpp"
#include "smoother.hpp"

using namespace std;

extern "C" {
int main_build(int argc, char *argv[]);
}

const string version = "v2.1.0";

int main(int argc, char **argv) {
  time_t t;
  time(&t);
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));

  // Print help if no arguments
  if (argc == 1) {
    cerr << MAIN_USAGE_MESSAGE << endl;
    // c->print_help("");
    exit(EXIT_FAILURE);
  }

  // Run ropebwt3 build if index mode
  if (strcmp(argv[1], "index") == 0) {
    spdlog::info(
        "We rely on 'ropebwt3 build' - just exposing it for convenience");
    main_build(argc - 1, argv + 1);
  } else {
    // Get configuration and run
    auto c = Configuration::getInstance();
    c->parse(argc, argv);
    if (c->verbose)
      spdlog::set_level(spdlog::level::debug);

    if (c->version) {
      cout << "SVDSS, " << version << endl;
      exit(EXIT_SUCCESS);
    }

    if (c->help) {
      c->print_help(argv[1]);
      exit(EXIT_SUCCESS);
    }

    if (strcmp(argv[1], "call") == 0) {
      if (c->reference == "" || c->bam == "" || c->sfs == "") {
        c->print_help(argv[1]);
        exit(EXIT_FAILURE);
      }
      auto caller = new Caller();
      caller->run();
    } else if (strcmp(argv[1], "search") == 0) {
      if (c->index == "" || (c->fastq == "" && c->bam == "")) {
        c->print_help(argv[1]);
        exit(EXIT_FAILURE);
      }
      auto pingpong = new PingPong();
      pingpong->search();
      // } else if (strcmp(argv[1], "assemble") == 0) {
      //   auto assembler = new Assembler();
      //   assembler->run();
    } else if (strcmp(argv[1], "smooth") == 0) {
      if (c->reference == "" || c->bam == "") {
        c->print_help(argv[1]);
        exit(EXIT_FAILURE);
      }
      Smoother().run();
    } else {
      c->print_help("");
      exit(EXIT_FAILURE);
    }
  }
  time_t s;
  time(&s);
  spdlog::info("All done! Runtime: {} seconds", s - t);

  return 0;
}
