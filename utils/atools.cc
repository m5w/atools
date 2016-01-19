#include <iostream>
#include <sstream>
#include <vector>

#include <queue>
#include <map>
#include <boost/program_options.hpp>
#include <getopt.h>
#include <unistd.h>
#include <memory>

#include "filelib.h"
#include "alignment_io.h"
#include "optional.h"

namespace po = boost::program_options;
using namespace std;

struct Command {
  virtual ~Command() {}
  virtual string Name() const = 0;

  // returns 1 for alignment grid output [default]
  // returns 2 if Summary() should be called [for AER, etc]
  virtual int Result() const { return 1; }

  virtual bool RequiresTwoOperands() const { return true; }
  virtual void Apply(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) = 0;
  void EnsureSize(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) {
    x->resize(max(a.width(), b.width()), max(a.height(), b.height()));
  }
  static bool Safe(const Array2D<bool>& a, int i, int j) {
    if (i >= 0 && j >= 0 && i < static_cast<int>(a.width()) && j < static_cast<int>(a.height()))
      return a(i,j);
    else
      return false;
  }
  virtual void Summary() { assert(!"Summary should have been overridden"); }
};

// compute fmeasure, second alignment is reference, first is hyp
struct FMeasureCommand : public Command {
  FMeasureCommand() : matches(), num_predicted(), num_in_ref() {}
  int Result() const { return 2; }
  string Name() const { return "fmeasure"; }
  bool RequiresTwoOperands() const { return true; }
  void Apply(const Array2D<bool>& hyp, const Array2D<bool>& ref, Array2D<bool>* x) {
    (void) x;   // AER just computes statistics, not an alignment
    unsigned i_len = ref.width();
    unsigned j_len = ref.height();
    for (unsigned i = 0; i < i_len; ++i) {
      for (unsigned j = 0; j < j_len; ++j) {
        if (ref(i,j)) {
          ++num_in_ref;
          if (Safe(hyp, i, j)) ++matches;
        } 
      }
    }
    for (unsigned i = 0; i < hyp.width(); ++i)
      for (unsigned j = 0; j < hyp.height(); ++j)
        if (hyp(i,j)) ++num_predicted;
  }
  void Summary() {
    if (num_predicted == 0 || num_in_ref == 0) {
      cerr << "Insufficient statistics to compute f-measure!\n";
      abort();
    }
    const double prec = static_cast<double>(matches) / num_predicted;
    const double rec = static_cast<double>(matches) / num_in_ref;
    cout << "P: " << prec << endl;
    cout << "R: " << rec << endl;
    const double f = (2.0 * prec * rec) / (rec + prec);
    cout << "F: " << f << endl;
  }
  int matches;
  int num_predicted;
  int num_in_ref;
};

struct DisplayCommand : public Command {
  string Name() const { return "display"; }
  bool RequiresTwoOperands() const { return false; }
  void Apply(const Array2D<bool>& in, const Array2D<bool>&, Array2D<bool>* x) {
    *x = in;
    cout << *x << endl;
  }
};

struct ConvertCommand : public Command {
  string Name() const { return "convert"; }
  bool RequiresTwoOperands() const { return false; }
  void Apply(const Array2D<bool>& in, const Array2D<bool>&, Array2D<bool>* x) {
    *x = in;
  }
};

struct InvertCommand : public Command {
  string Name() const { return "invert"; }
  bool RequiresTwoOperands() const { return false; }
  void Apply(const Array2D<bool>& in, const Array2D<bool>&, Array2D<bool>* x) {
    Array2D<bool>& res = *x;
    res.resize(in.height(), in.width());
    for (unsigned i = 0; i < in.height(); ++i)
      for (unsigned j = 0; j < in.width(); ++j)
        res(i, j) = in(j, i);
  }
};

struct IntersectCommand : public Command {
  string Name() const { return "intersect"; }
  bool RequiresTwoOperands() const { return true; }
  void Apply(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) {
    EnsureSize(a, b, x);
    Array2D<bool>& res = *x;
    for (unsigned i = 0; i < a.width(); ++i)
      for (unsigned j = 0; j < a.height(); ++j)
        res(i, j) = Safe(a, i, j) && Safe(b, i, j);
  }
};

struct UnionCommand : public Command {
  string Name() const { return "union"; }
  bool RequiresTwoOperands() const { return true; }
  void Apply(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) {
    EnsureSize(a, b, x);
    Array2D<bool>& res = *x;
    for (unsigned i = 0; i < res.width(); ++i)
      for (unsigned j = 0; j < res.height(); ++j)
        res(i, j) = Safe(a, i, j) || Safe(b, i, j);
  }
};

struct RefineCommand : public Command {
  RefineCommand() {
    neighbors_.push_back(make_pair(1,0));
    neighbors_.push_back(make_pair(-1,0));
    neighbors_.push_back(make_pair(0,1));
    neighbors_.push_back(make_pair(0,-1));
  }
  bool RequiresTwoOperands() const { return true; }

  void Align(unsigned i, unsigned j) {
    res_(i, j) = true;
    is_i_aligned_[i] = true;
    is_j_aligned_[j] = true;
  }

  bool IsNeighborAligned(int i, int j) const {
    for (unsigned k = 0; k < neighbors_.size(); ++k) {
      const int di = neighbors_[k].first;
      const int dj = neighbors_[k].second;
      if (Safe(res_, i + di, j + dj))
        return true;
    }
    return false;
  }

  bool IsNeitherAligned(int i, int j) const {
    return !(is_i_aligned_[i] || is_j_aligned_[j]);
  }

  bool IsOneOrBothUnaligned(int i, int j) const {
    return !(is_i_aligned_[i] && is_j_aligned_[j]);
  }

  bool KoehnAligned(int i, int j) const {
    return IsOneOrBothUnaligned(i, j) && IsNeighborAligned(i, j);
  }

  typedef bool (RefineCommand::*Predicate)(int i, int j) const;

 protected:
  void InitRefine(
      const Array2D<bool>& a,
      const Array2D<bool>& b) {
    res_.clear();
    EnsureSize(a, b, &res_);
    in_.clear(); un_.clear(); is_i_aligned_.clear(); is_j_aligned_.clear();
    EnsureSize(a, b, &in_);
    EnsureSize(a, b, &un_);
    is_i_aligned_.resize(res_.width(), false);
    is_j_aligned_.resize(res_.height(), false);
    for (unsigned i = 0; i < in_.width(); ++i)
      for (unsigned j = 0; j < in_.height(); ++j) {
        un_(i, j) = Safe(a, i, j) || Safe(b, i, j);
        in_(i, j) = Safe(a, i, j) && Safe(b, i, j);
        if (in_(i, j)) Align(i, j);
    }
  }
  // "grow" the resulting alignment using the points in adds
  // if they match the constraints determined by pred
  void Grow(Predicate pred, bool idempotent, const Array2D<bool>& adds) {
    if (idempotent) {
      for (unsigned i = 0; i < adds.width(); ++i)
        for (unsigned j = 0; j < adds.height(); ++j) {
          if (adds(i, j) && !res_(i, j) &&
              (this->*pred)(i, j)) Align(i, j);
        }
      return;
    }
    set<pair<int, int> > p;
    for (unsigned i = 0; i < adds.width(); ++i)
      for (unsigned j = 0; j < adds.height(); ++j)
        if (adds(i, j) && !res_(i, j))
          p.insert(make_pair(i, j));
    bool keep_going = !p.empty();
    while (keep_going) {
      keep_going = false;
      set<pair<int, int> > added;
      for (set<pair<int, int> >::iterator pi = p.begin(); pi != p.end(); ++pi) {
        if ((this->*pred)(pi->first, pi->second)) {
          Align(pi->first, pi->second);
          added.insert(make_pair(pi->first, pi->second));
          keep_going = true;
        }
      }
      for (set<pair<int, int> >::iterator ai = added.begin(); ai != added.end(); ++ai)
        p.erase(*ai);
    }
  }
  Array2D<bool> res_;  // refined alignment
  Array2D<bool> in_;   // intersection alignment
  Array2D<bool> un_;   // union alignment
  vector<bool> is_i_aligned_;
  vector<bool> is_j_aligned_;
  vector<pair<int,int> > neighbors_;
};

struct DiagCommand : public RefineCommand {
  DiagCommand() {
    neighbors_.push_back(make_pair(1,1));
    neighbors_.push_back(make_pair(-1,1));
    neighbors_.push_back(make_pair(1,-1));
    neighbors_.push_back(make_pair(-1,-1));
  }
};

struct GDCommand : public DiagCommand {
  string Name() const { return "grow-diag"; }
  void Apply(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) {
    InitRefine(a, b);
    Grow(&RefineCommand::KoehnAligned, false, un_);
    *x = res_;
  }
};

struct GDFCommand : public DiagCommand {
  string Name() const { return "grow-diag-final"; }
  void Apply(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) {
    InitRefine(a, b);
    Grow(&RefineCommand::KoehnAligned, false, un_);
    Grow(&RefineCommand::IsOneOrBothUnaligned, true, a);
    Grow(&RefineCommand::IsOneOrBothUnaligned, true, b);
    *x = res_;
  }
};

struct GDFACommand : public DiagCommand {
  string Name() const { return "grow-diag-final-and"; }
  void Apply(const Array2D<bool>& a, const Array2D<bool>& b, Array2D<bool>* x) {
    InitRefine(a, b);
    Grow(&RefineCommand::KoehnAligned, false, un_);
    Grow(&RefineCommand::IsNeitherAligned, true, a);
    Grow(&RefineCommand::IsNeitherAligned, true, b);
    *x = res_;
  }
};

map<string, shared_ptr<Command> > commands;

class Help {
public:
  static const Help &help();
  void operator()() const;
private:
  Help();
  std::ostringstream dcmdline_options;
};

const Help &Help::help() { static Help TheHelp; return TheHelp; }

void Help::operator()() const {
  std::cout << dcmdline_options.str();
}

Help::Help() : dcmdline_options() {
  dcmdline_options <<
"\n"
"Configuration options:\n"
"  -i [ --input_1 ] arg            [REQUIRED] Alignment 1 file, - for STDIN\n"
"  -j [ --input_2 ] arg            Alignment 2 file, - for STDIN\n"
"  -c [ --command ] arg (=convert) Operation to perform:";

  for (auto Command_ : commands)
    dcmdline_options << ' ' << Command_.first;
}

struct option longopts[] = {
  {"command", required_argument, 0, 'c'},
  {"help", no_argument, 0, 'h'},
  {"input_1", required_argument, 0, 'i'},
  {"input_2", required_argument, 0, 'j'}
};

int indexptr;

void set_indexptr(const int val) {
  if (val == longopts[indexptr].val)
    return;

  for (std::size_t longopts_Index = 0;; ++longopts_Index) {
    const int val_ = longopts[longopts_Index].val;

    if (val_ == 0)
      break;

    if (val_ == val) {
      indexptr = longopts_Index;
      return;
    }
  }
}

std::string option_string(const struct option option_) {
  std::ostringstream option_string_;
  option_string_ << '-' << static_cast<char>(option_.val) << ", --"
                 << option_.name;
  return option_string_.str();
}

std::string option_string(const int indexptr_) {
  return option_string(longopts[indexptr_]);
}

std::string option_string() { return option_string(indexptr);
}

class UnexpectedOption : public std::exception {
public:
  UnexpectedOption(const std::ostringstream &what_);
  ~UnexpectedOption() throw();
  const char *what() const throw();
private:
  const std::string what_;
};

UnexpectedOption::UnexpectedOption(const std::ostringstream &what_)
    : what_(what_.str()) {}

UnexpectedOption::~UnexpectedOption() throw() {}

const char *UnexpectedOption::what() const throw() { return what_.c_str(); }

void optionCase(Apertium::Optional<std::string> &optarg_) {
  if (optarg_) {
    std::ostringstream what_;
    what_ << "unexpected option " << option_string() << " following option "
          << option_string();
    throw UnexpectedOption(what_);
  }

  optarg_ = std::string(optarg);
}

Apertium::Optional<std::string> command, input_1, input_2;

class getopt_long_Exception {};

void InitCommandLine(unsigned argc, char** argv) {

  while (true) {
    int getopt_long_ = getopt_long(argc, argv, "c:hi:j:", longopts, &indexptr);

    if (getopt_long_ == -1)
      break;

    set_indexptr(getopt_long_);

    switch (getopt_long_) {
    case 'c':
      optionCase(command);
      break;
    case 'h':
      Help::help()();
      exit(0);
    case 'i':
      optionCase(input_1);
      break;
    case 'j':
      optionCase(input_2);
      break;
    default:
      Help::help()();
      throw getopt_long_Exception();
    }
  }

  if (!command)
    command = std::string("convert");

  if (!input_1) {
    Help::help()();
    exit(1);
  }

  if (commands.count(*command) == 0) {
    cerr << "Don't understand command: " << *command << endl;
    exit(1);
  }
  if (commands[*command]->RequiresTwoOperands()) {
    if (!input_2) {
      cerr << "Command '" << *command << "' requires two alignment files\n";
      exit(1);
    }
    if (*input_1 == "-" && *input_2 == "-") {
      cerr << "Both inputs STDIN, reading PAIRS of lines\n";
    }
  } else {
    if (input_2) {
      cerr << "Command '" << *command << "' requires only one alignment file\n";
      exit(1);
    }
  }
}

template<class C> static void AddCommand() {
  C* c = new C;
  commands[c->Name()].reset(c);
}

int main(int argc, char **argv) {
  AddCommand<ConvertCommand>();
  AddCommand<DisplayCommand>();
  AddCommand<InvertCommand>();
  AddCommand<IntersectCommand>();
  AddCommand<UnionCommand>();
  AddCommand<GDCommand>();
  AddCommand<GDFCommand>();
  AddCommand<GDFACommand>();
  AddCommand<FMeasureCommand>();
  InitCommandLine(argc, argv);
  Command& cmd = *commands[*command];
  shared_ptr<ReadFile> rf1(new ReadFile(*input_1));
  shared_ptr<ReadFile> rf2;
  if (cmd.RequiresTwoOperands())
    rf2.reset(new ReadFile(*input_2));
  istream* in1 = rf1->stream();
  istream* in2 = NULL;
  if (rf2) in2 = rf2->stream();
  while(*in1) {
    string line1;
    string line2;
    getline(*in1, line1);
    if (in2) {
      getline(*in2, line2);
      if ((*in1 && !*in2) || (*in2 && !*in1)) {
        cerr << "Mismatched number of lines!\n";
        exit(1);
      }
    }
    if (line1.empty() && !*in1) break;
    shared_ptr<Array2D<bool> > out(new Array2D<bool>);
    shared_ptr<Array2D<bool> > a1 = AlignmentIO::ReadPharaohAlignmentGrid(line1);
    if (in2) {
      shared_ptr<Array2D<bool> > a2 = AlignmentIO::ReadPharaohAlignmentGrid(line2);
      cmd.Apply(*a1, *a2, out.get());
    } else {
      Array2D<bool> dummy;
      cmd.Apply(*a1, dummy, out.get());
    }
    
    if (cmd.Result() == 1) {
      AlignmentIO::SerializePharaohFormat(*out, &cout);
    }
  }
  if (cmd.Result() == 2)
    cmd.Summary();
  return 0;
}

